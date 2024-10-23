function  [Nd,El,OHC,IHC,MP,FLD] = initiate(MP)

    global coef
    opt_fin = 0;  % 0: generate new fluid-domain mesh, 1: use existing mesh grid
    opt_MP_set = 1;         % 0: the MP set before 2021.04
    
    dim = 3;    % dimension 2 or 3-D

    opt_lgradient = 1; % longitudinal gradient
    if MP.vMC == 1, opt_lgradient = 1; end % to control the gradient for testing purposes
        
    [Nd, El, MP] = model_3D(MP, dim, opt_lgradient, opt_MP_set);
    
    MP.opt_lgradient = opt_lgradient;

    if MP.vMC == 1
        MP.EP = 0; % endocochlear potential in mV
    else
        MP.EP = 90;
    end
    
    MP.EK = 75; % cell equilibrium potential in mV,note the sign, MSA: is this even used?!

    Nd = set_boundary_conditions(MP, Nd, opt_MP_set);

    [OHC, IHC] = model_OHC_2state_circuit(MP,El);

    FLD = model_FLD_cochlea(MP,MP.H,Nd);
    if opt_fin == 1
        input = load('.\VC_code\hinput\mesh_12mm_cochlea_quad_4.mat','FLD');
        if isfield(input.FLD,'p')
            FLD.h0 = input.FLD.h0;
            FLD.HOC = input.FLD.HOC;
            FLD.L = input.FLD.L;
            FLD.Nd.x = input.FLD.p;
            FLD.El.node = input.FLD.t;
            FLD.El.N = size(FLD.El.node,1);
            FLD.Nd.N = size(FLD.Nd.x,1);
            FLD.El.type = 3;
        end                
    else          
        FLD = FLD_mesh(MP,MP.H,FLD);
        % the following lines for when quadrilateral elements are used
        if eq(FLD.El.type,4) 
            FLD.El.NQ = 4;
%             set_GaussQuad_as_a_global_variable(FLD.El.NQ);
            [xi, eta, ~] = set_Gauss_local_variables(2,2,FLD.El.NQ);
            FLD.El = evaluate_interp_func(FLD.Nd,FLD.El,xi,eta);
        end
    end
    
    if MP.vMC == 1
        FLD = set_FEfluid_BCs_MC(MP,FLD);
    else
        FLD = set_FEfluid_BCs(MP,FLD);
    end

    if MP.opt_cfld == 1
        CFLD = model_CFLD(MP,Nd);
        CFLD = CFLD_mesh_2D(MP,Nd,CFLD);
        
        CFLD.El.NQ = 4;
        CFLD = set_CFLD_BC_2D(CFLD);
        [xi, eta, ~] = set_Gauss_local_variables(2,2,CFLD.El.NQ);
        CFLD.El = evaluate_interp_func(CFLD.Nd,CFLD.El,xi,eta);

        FLD.CFLD = CFLD;
    end

    if MP.opt_cfld == 1      
        % permeability condition
        alpha = MP.p_alpha;
        kk = alpha;
        MP.kk = kk;

    end

    % Set damping
    alph_mltplr = 1;
    ac_a = 3*0.06*alph_mltplr;      % testing
    ac_b = 2*0.001*alph_mltplr;     % testing
    beta_multiplier = coef.beta_multiplier;
    beta_ca = 10*beta_multiplier;
    beta_cb = 20*beta_multiplier;
    gradb = log(beta_ca/beta_cb)/8000;
    grad = log(ac_a/ac_b)/8000;
    MP.alpha_cgrad = grad;      % longitudinal gradient of damping
    MP.beta_cgrad = gradb;
    if MP.loc >6
        MP.beta_c = beta_ca*exp(gradb*(MP.loc-10)*1e3);
        MP.alpha_c = ac_a*exp(grad*(MP.loc-10)*1e3);    % damping factor, C = alpha*K
    elseif MP.loc <= 6
        MP.alpha_c = ac_b*exp(grad*(MP.loc-2)*1e3);     % damping factor, C = alpha*K
        MP.beta_c = beta_cb*exp(gradb*(MP.loc-2)*1e3);
    end

    global coef
    if ~isempty(coef)
        if isfield(coef,'BM_damping_coef')
            MP.BM_damping_coef = coef.BM_damping_coef;
        else
            MP.BM_damping_coef = beta_multiplier;
        end
        if isfield(coef,'RL_damping_coef')
            MP.RL_damping_coef = coef.RL_damping_coef;
        else
            MP.RL_damping_coef = beta_multiplier;
        end
        if isfield(coef,'TM_damping_coef')
            MP.TM_damping_coef = coef.TM_damping_coef;
        else
            MP.TM_damping_coef = beta_multiplier;
        end
        if isfield(coef,'DC_damping_coef')
            MP.DC_damping_coef = coef.DC_damping_coef;
        else
            MP.DC_damping_coef = beta_multiplier;
        end        
        if isfield(coef,'PC_damping_coef')
            MP.PC_damping_coef = coef.PC_damping_coef;
        else
            MP.PC_damping_coef = beta_multiplier;
        end        
    else
        MP.BM_damping_coef = beta_multiplier;
        MP.RL_damping_coef = beta_multiplier;
        MP.TM_damping_coef = beta_multiplier;
        MP.DC_damping_coef = beta_multiplier;
        MP.PC_damping_coef = beta_multiplier;
    end    

    % initiating degree of freedom structure
    
    MP.dof(1).name = 'scala pressure';
    MP.dof(1).n = FLD.tdof;
    MP.dof(1).pdof = 1:MP.dof(1).n;
    MP.dof(1).dof_r = 1:FLD.rdof;
    last_dof_r = MP.dof(1).dof_r(end);

    MP.dof(2).name = 'structural motion';
    MP.dof(2).n = Nd.tdof;
    MP.dof(2).udof = MP.dof(1).n + (1:MP.dof(2).n);
    MP.dof(2).dof_r = last_dof_r + (1:Nd.rdof);
    last_dof_r = MP.dof(2).dof_r(end);

    nz = length(MP.xx);
    if MP.ind_MET == 2
        nv = 5; % 5 for node voltage
        MP.dof(3).name = '5 node circuit';
        MP.dof(3).n = nv*nz;
        MP.dof(3).edof = sum([MP.dof([1,2]).n]) + (1:MP.dof(3).n);
        MP.dof(3).dof_r = last_dof_r + (1:MP.dof(3).n);
        last_dof_r = MP.dof(3).dof_r(end);

        MP.dof(4).name = 'adaptation x_a';
        MP.dof(4).n = nz;
        MP.dof(4).adof = sum([MP.dof((1:3)).n]) + (1:MP.dof(4).n);
        MP.dof(4).dof_r = last_dof_r + (1:MP.dof(4).n);
        last_dof_r = MP.dof(4).dof_r(end);

        MP.dof(5).name = 'channel open probability';
        MP.dof(5).n = nz;
        MP.dof(5).odof = sum([MP.dof((1:4)).n]) + (1:MP.dof(5).n);
        MP.dof(5).dof_r = last_dof_r + (1:MP.dof(5).n);
        last_dof_r = MP.dof(5).dof_r(end);
    else
        error('10 state channel is not implemented in this code');
    end
    
    if MP.opt_cfld == 1
        MP.dof(6).name = 'Corti fluid';
        MP.dof(6).n = CFLD.tdof;
        MP.dof(6).cdof = sum([MP.dof((1:5)).n]) + (1:MP.dof(6).n);
        MP.dof(6).dof_r = last_dof_r + (1:CFLD.rdof);
        last_dof_r = MP.dof(6).dof_r(end);

        MP.dof(7).name = "Peristalsis";
        MP.dof(7).n = CFLD.Nd.nz;
        MP.dof(7).ldof = sum([MP.dof((1:6)).n]) + (1:MP.dof(7).n);
        BC = reshape(FLD.CFLD.BC,4,FLD.CFLD.Nd.N);
        p_BC = BC(2,FLD.CFLD.Nd.ind_Radi);
        ldof_r = sum(BC(2,FLD.CFLD.Nd.ind_Radi));
        MP.dof(7).dof_r = last_dof_r + (1:ldof_r);
        last_dof_r = MP.dof(7).dof_r(end);
    else
        MP.dof(6).n = 0;
        MP.dof(7).n = 0;
    end

    MP.dof(8).name = 'Lagrange constraint';
    MP.dof(8).n = 3*0.5*numel(Nd.Master_Slave); % [dx,dy,dz] only
    MP.dof(8).ldof = sum([MP.dof((1:7)).n]) + (1:MP.dof(8).n);
    nd = Nd.Master_Slave(:,1);
    rdof = (nd-1)*6 + (1:3);
    rdof = Nd.mapF2R(rdof);
    l_BC = logical(reshape(transpose(rdof),MP.dof(8).n,1));
    ldof_r = sum(l_BC);
    MP.dof(8).dof_r = last_dof_r + (1:ldof_r);   

    if MP.opt_cfld == 1
        MP.BC = [FLD.BC; Nd.BC; true(MP.dof(3).n,1); true(MP.dof(4).n,1);...
        true(MP.dof(5).n,1); CFLD.BC; p_BC(:); l_BC];
        MP.BC_psv = [FLD.BC; Nd.BC; false(MP.dof(3).n,1); false(MP.dof(4).n,1);...
        false(MP.dof(5).n,1); CFLD.BC; p_BC(:); l_BC];
    else
        MP.BC = [FLD.BC; Nd.BC; true(MP.dof(3).n,1); true(MP.dof(4).n,1);...
        true(MP.dof(5).n,1); l_BC];
        MP.BC_psv = [FLD.BC; Nd.BC; false(MP.dof(3).n,1); false(MP.dof(4).n,1);...
        false(MP.dof(5).n,1); l_BC];
    end

    MP.tdof = length(MP.BC);
    MP.rdof = sum(MP.BC);

    % for reducing matrix and rhs from reduced active to reduced passive
    MP.act2psv = true(MP.rdof,1);
    MP.act2psv([MP.dof(3).dof_r,MP.dof(4).dof_r,MP.dof(5).dof_r]) = false;
    MP.rdof_psv = sum(MP.act2psv);
end % of function initiate()
