function [Uf,AMat] = NLFD_model(MP,fMP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0,Uf,opt_lin)

    tol = 1e-2;

    HB = [OHC.HB];
    po0 = transpose([HB.po0]);

    if opt_lin
        % po wrt xHB and xa: po - d(po)dxHB*xHB - d(po)dxa*xa = po_rest - d(po)dxHB*xHB_rest - d(po)dxa*xa_rest
        AMat = make_lin_Aou_Aoa(AMat, MP, Nd, El, OHC); % matrices for d(po)dxHB and d(po)dxa
        % stereocilia membrane conductance (Gs) wrt po: after taylor expansion on Gs_max*po*V
        % iwCe*V + Ge*V + [Gs_max*V0*p -> (Aeo)] + [Gs_max*po0*V -> included in Ge] = I + [Gs_max*po0*V0 -> included in loading condition] 
        AMat.Aeo = make_lin_Aeo(MP, OHC);
        V_lin = AMat.Aeo*po0;

        % append A0
        udof = MP.dof(2).dof_r;
        edof = MP.dof(3).dof_r;
        adof = MP.dof(4).dof_r;
        odof = MP.dof(5).dof_r;

        A0(odof,udof) = A0(odof,udof) + AMat.Aou;
        A0(odof,adof) = A0(odof,adof) + AMat.Aoa;
        A0(edof,odof) = A0(edof,odof) + AMat.Aeo;
    else
        V_lin = zeros(MP.dof(3).n,1);
    end

    [Ff, Uf, FLD] = loading_condition(MP,fMP,FLD,V_lin,Ie,po0,Uf); % pressure and Ie applied
    % Rf is a reduced dof loading vector
    % Uf is updated to contain input stapes pressure value
    
    icount = 0;

    dA = cell(1,fMP.L); % crude way to initiate the decomposition cell array of LHS

    if ~opt_lin
        if (MP.key_fOHC || MP.key_fMET)
            h = initiate_plots;       
        end
        fprintf(1,'\n  harmonics iterations relative error:\n'); % initiate error report
    end

    while true
       
        icount = icount + 1;
        dtau = MP.dtau;

        Uf_n = Uf;
        if (MP.key_fOHC || MP.key_fMET) && ~opt_lin
            harmonics_threshold = 20;
        else
            harmonics_threshold = 0;
        end

        [fMP, Uf, dA] = evaluate_Uf(dtau,MP,fMP,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,V0,po0,Uf,dA,Ff,icount,harmonics_threshold,opt_lin);
      
        if (MP.key_fOHC || MP.key_fMET) && ~opt_lin % active cochlea condition
            plot_output(h,MP,fMP,Nd,El,AMat,Uf,po0,V0)
        drawnow
        end

        if ~(MP.key_fOHC || MP.key_fMET || MP.opt_cfld) || opt_lin
            if (MP.key_fOHC || MP.key_fMET) && ~MP.vMC
                Uf = correct_MET_saturation(dtau,MP,fMP,Nd,El,OHC,FLD,V_lin,A0,AMat,Ce,Ge,Ie,V0,po0,Uf,dA,icount,opt_lin);
            end
            break;
        end

        if (MP.key_fOHC || MP.key_fMET)
            g_res_val = evaluate_err(fMP,MP,Nd,PoI,Uf_n,Uf,icount,h);
        else
            g_res_val = evaluate_err(fMP,MP,Nd,PoI,Uf_n,Uf,icount);
        end
%         if (g_res_val < tol) || (icount == 10) || (g_res_val > 1e13)
        if ((g_res_val < tol) || (g_res_val > 1e13))
            break;
        end
    end
end

function Uf = correct_MET_saturation(dtau,MP,fMP,Nd,El,OHC,FLD,V_lin,A0,AMat,Ce,Ge,Ie,V0,po0,Uf,dA,icount,opt_lin)

    odof = MP.dof(5).odof;
    po = Uf(odof,fMP.ind); 
    % assuming po_rest is the same for all OHCs
    po_p = zeros(MP.dof(5).n,fMP.L); % prescribed saturated open channel probability
    col = fMP.ind(2);
    z_col = fMP.ind(1);
    % find saturated channels
    factor = 0.8;    
    if median(po(:,1)) < 0.5        
        idx = (abs(po(:,2)) > factor*po(:,1));
        po_p(idx,col) = po(idx,2);
        idx_sat = (abs(po(:,2)) > po(:,1));
        po_p(idx_sat,col) = po(idx_sat,2)./abs(po(idx_sat,2)).*po(idx_sat,1);
        po_p(idx,col) = po_p(idx,col)./abs(po_p(idx,col)).*smooth(abs(po_p(idx,col)),0.25);
        po_p(:,z_col) = po(:,1);
    else
        idx = ((abs(po(:,2)) + po(:,1)) > factor);
        po_p(idx,col) = po(idx,2)./abs(po(idx,2)).*(1-po(idx,1));
        po_p(:,z_col) = po(:,1);
    end

    %update the po prescribed vector for complex values
    po_p(:,fMP.L-col+2) = conj(po_p(:,col));

    % update the output for prescribed values
    Uf(odof,:) = po_p;

    % reduce A0 to account for prescribed po! this is required since the code was designed to prematurely implement matrix reduction
    rdof = true(MP.rdof,1);
    rdof(MP.dof(5).dof_r(idx)) = false;
    A0 = A0(rdof,rdof);

    % update the unknown dof vector
    MP.BC(odof,1) = ~idx;
    MP.dof(5).dof_r = MP.dof(4).dof_r(end) + (1:sum(~idx)); % odof rdof updated
    
    if MP.opt_cfld == 1
        MP.dof(6).dof_r = MP.dof(5).dof_r(end) + (1:FLD.CFLD.rdof);
        MP.dof(7).dof_r = MP.dof(6).dof_r(end) + (1:MP.dof(7).n);
    end
    MP.rdof = sum(MP.BC);

    % reproduce the loading condition vector, this time w input SPL instead of zero
    [Ff, Uf, FLD] = loading_condition(MP,fMP,FLD,V_lin,Ie,po0,Uf); % pressure and Ie applied

    % update the RHS for prescribed po via Auo and Aeo matrices
    row = MP.dof(2).dof_r;
    Ff(row,:) = Ff(row,:) - AMat.Auo*po_p;
    row = MP.dof(3).dof_r;
    Ff(row,:) = Ff(row,:) - AMat.Aeo*po_p;
    
    [~, Uf, ~] = evaluate_Uf(dtau,MP,fMP,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,V0,po0,Uf,dA,Ff,icount,[],opt_lin);

end

function [fMP, Uf, dA] = evaluate_Uf(dtau,MP,fMP,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,V0,po0,Uf,dA,Ff,icount,threshold,opt_lin)

    MP.dtau = dtau;
    key_fOHC = MP.key_fOHC;
    key_fMET = MP.key_fMET;
    opt_cfld = MP.opt_cfld;
    R = Ff;
    clear Ff;

    %%%%%% gather the data needed for evaluating the nonlinear terms
    if (opt_cfld == 1) && ~opt_lin % Corti fluid nonlinear terms
        CFLD = FLD.CFLD;
        convection = zeros(CFLD.El.N,2*CFLD.El.NQ, fMP.L);        
        % extract required terms for Corti fluid nonlinearity
        idx = reshape(MP.dof(6).cdof,4,[]);    
        CUf = Uf(idx(1,:)',:);       
        CVf = Uf(idx(2,:)',:);        
        CU = inverse_fourier_transform(CUf);
        CV = inverse_fourier_transform(CVf);
        % Corti flux is updated from last iteration since it is 
        % being used in the rhs vector of the pseudo time iterative method
        MP.bl_th = sqrt(MP.nu/fMP.freq);
        for ti = 1:fMP.L
            if opt_cfld == 1
                convection(:,:,ti) = ...
                    evaluate_CFLD_nonlinear_terms(FLD, CU(:,ti), ...
                    CV(:,ti));
            end      
        end
        convection = fourier_transform(convection);
        if icount > 1
            fMP = identify_harmonics(fMP,convection,20);
        end
        CFLD.convection = convection;
        FLD.CFLD = CFLD; 
    end

    if (key_fOHC || key_fMET) && ~opt_lin % active cochlea condition
        GV_nonlin = zeros(MP.dof(3).n, fMP.L); % initiate the nonlinear component of G*V in the circuit model            
        % extract outer hair bundle deflection from element ANK
        xxHB = HB_deflection(MP,El,Uf,60,icount);
        idx = MP.dof(4).adof;
        xa = inverse_fourier_transform(Uf(idx,:));
        idx = MP.dof(5).odof;        
        po_f = Uf(idx,:);
        po = inverse_fourier_transform(po_f);
        % extract voltage vector
        idx = MP.dof(3).edof;        
        voltage_f = Uf(idx,:);
        voltage = inverse_fourier_transform(Uf(idx,:));
        HB = [OHC.HB];
        po_rhs = zeros(MP.nz, fMP.L);        
        parfor ti = 1:fMP.L
            po_rhs(:,ti) = channel_2state(HB, xxHB(:,ti), xa(:,ti));                            
        end        
        po_f_rhs = fourier_transform(po_rhs);
        idx = (fMP.ffsym == 0);
        po0 = po_f(:,idx);
        parfor ti = 1:fMP.L
            GV_nonlin(:,ti) = evaluate_GV_nonlin(MP,OHC,voltage(:,ti),po(:,ti),po0);
        end
        GV_nonlin_f = fourier_transform(GV_nonlin);

        if icount > 1
            nonlin = po_f_rhs;
            nonlin(:,fMP.ffsym == 0) = 0;
            fMP = identify_harmonics(fMP,nonlin,threshold);
        end
        clear 'HB' 'xxHB' 'xa' 'voltage' 'po' 'po_rhs' 'GV_nonlin';
    end

    % the reason forces are calculated mid iteration is that they might
    % change due to changes in resting po and V, might be able to skip this
    % if resting po was 0.5 from the start
    if key_fOHC || key_fMET % active cochlea condition
        Auo = AMat.Auo; % used to evaluate fMET rhs at rest
        Aue = AMat.Aue; % used to evaluate fOHC rhs at rest
    end

    idx = MP.dof(2).udof;
    uu_f = Uf(idx,fMP.ind);       
    idx = MP.dof(1).pdof;
    p_f = Uf(idx,fMP.ind);

    R = R(:,fMP.ind);

    for i = 1:length(fMP.ind) % loop over stim frequencies and their harmonics

        fi = fMP.ind(i);
        wk = 2*pi*fMP.ffsym(fi);       
        
        % RHS for CFLD nonlinearity is evaluated along with CLFD matrices
        if opt_cfld == 1
            FLD = make_CFLD_mat_nl(fi, MP, FLD, Uf(MP.dof(6).cdof,fi), wk,opt_lin);
            b = zeros(MP.tdof,1);
            idx = MP.dof(6).cdof;
            b(idx) = FLD.CFLD.b;
            R(:,i) = R(:,i) + b(MP.BC,1);
            AMat.L = make_peristalsis_effect(MP, FLD, wk);
        end

        [A1, AMat] = assemble_A1(Nd, FLD, MP, AMat, wk, Ce, Ge, opt_lin);
        if (icount == 1) || isempty(dA{fi}) % condition to make the new matrices for the rising harmonics
            A = A0 + A1;
            if ~( key_fOHC || key_fMET ) % active cochlea condition
                A = A(MP.act2psv,MP.act2psv);
            end
            dA{fi} = decompose_matrix(A);
        end  

        if key_fOHC || key_fMET % active cochlea condition
            if ~opt_lin
                % RHS for mechanotransduction nonlinearity
                b = make_MET_RHS(opt_lin,AMat,wk,MP,V0,po0,Auo,Aue,po_f(:,fi),...
                    voltage_f(:,fi),po_f_rhs(:,fi),GV_nonlin_f(:,fi));
            else
                b = make_MET_RHS(opt_lin,AMat,wk,MP,V0,po0,Auo,Aue);
            end
            R(:,i) = R(:,i) + b;
        end

        if MP.opt_structure_relaxation
            b = make_structure_RHS(MP,Nd,uu_f(:,i));
            R(:,i) = R(:,i) + b;
        end

        if MP.opt_fluid_relaxation
            b = make_scalae_RHS(MP,FLD,p_f(:,i));
            R(:,i) = R(:,i) + b;
        end
    end % end of assembly
    
    % prepare for Ax = b
    if key_fOHC || key_fMET % active cochlea condition
        row = MP.rdof;
    else
        row = MP.rdof_psv;
        R = R(MP.act2psv,:);
    end
    a = zeros(row,length(fMP.ind),'like',1j);

    % Ax = b
   
    for i = 1:length(fMP.ind) % loop over stim frequencies and their harmonics
        fi = fMP.ind(i);
        a(:,i) = dA{fi}\R(:,i);
    end

    % update Uf and relaxation
    ind = false(1,fMP.L);
    ind(fMP.ind) = true;
    if opt_lin 
        ind(fMP.ind(1)) = false;
        a(:,1) = [];
    end % do not update the DC on linear solution
    if key_fOHC || key_fMET % active cochlea condition
        Uf(MP.BC,ind) = a;
    else
        Uf(MP.BC_psv,ind) = a;
    end
    ind_conj = false(1,fMP.L);
    ind_conj(fMP.L-fMP.ind+2) = true;
    if opt_lin
        ind_conj(fMP.L-fMP.ind(1)+2) = false;
    end
    Uf(:,ind_conj) = fliplr(conj(Uf(:,ind)));    

end
