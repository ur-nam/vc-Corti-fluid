function [A1, G] = assemble_A1(Nd, FLD, MP, G, ww, Ce, Ge, opt_lin)

if MP.Visc == 0
    [A1, G] = assemble_inviscid(Nd, FLD, MP, G, ww, Ce, Ge, opt_lin);
end

end

function [A1, G] = assemble_inviscid(Nd, FLD, MP, G, ww, Ce, Ge, opt_lin)

% % function A1 = assemble_A1(Nd, El, OHC, MP, G, ww)
% %
% % Assembe frequency-dependent submatrices that must be updated with
% % stimulating frequency
% %
    
    if ~opt_lin
        dtau_inv = 1/MP.dtau;
    else
        dtau_inv = 0;
        
    end

    opt_cfld = MP.opt_cfld;
    tdof = MP.tdof;
    rdof = MP.rdof;

    pdof = MP.dof(1).dof_r;
    udof = MP.dof(2).dof_r;
    edof = MP.dof(3).dof_r;
    adof = MP.dof(4).dof_r;
    odof = MP.dof(5).dof_r;
    
    mdensity = 1.6e-4; % mdensity = nzmax(C1)/tdof^2, fraction of non-zero element, This value depend on FE model's connectivity
    N_max_nonzero = round(mdensity*tdof^2); % max number of non-zero elements.  Try N_max_nonzero = nzmax(C)
    
    A1 = spalloc(rdof,rdof,N_max_nonzero);
    

    Apu = (1i*ww)^2*G.Apu_O;
    A1(pdof, udof) = Apu;

    Auu = (1i*ww)^2*G.M + (1i*ww).*G.C + G.K;
    if MP.opt_structure_relaxation
        Auu = dtau_inv*speye(Nd.rdof) + Auu;
    end

    A1(udof, udof) = Auu;

    if MP.opt_fluid_relaxation
        A1(pdof, pdof) = dtau_inv*speye(FLD.rdof);
    end

    if MP.key_fOHC || MP.key_fMET % active cochlea condition

%         Aee = dtau_inv*speye(MP.dof(3).n) + 1i*ww*Ce + Ge;
        Aee = 1i*ww*Ce + Ge;
        G.Aee = Aee;
        A1(edof, edof) = Aee*(dtau_inv + 1);
    
        Aaa = 1i*ww*G.Iaa;
        A1(adof, adof) = Aaa;
        
        if ~isempty(odof)
            Aoo = speye(length(odof))*(1 + dtau_inv);
            A1(odof, odof) = Aoo;
        end
    end

%     Qpp = G.Qpp_O*exp(-MP.p_beta*ww);
%     A1(pdof, pdof) = A1(pdof, pdof) - (1i*ww)*Qpp;
    
    if opt_cfld == 1
        CFLD = FLD.CFLD;
        cdof = MP.dof(6).dof_r;
        Acc = CFLD.Acc;
        A1(cdof, cdof) = A1(cdof, cdof) + Acc;
%         Qpc = G.Qpc_O*exp(-MP.p_beta*ww);
%         A1(pdof, cdof) = A1(pdof, cdof) + (1i*ww)*Qpc;

        % Lagrange multiplier formulation for area change to radial motion
        % of Corti coupling
        ldof = MP.dof(7).dof_r;
        dof = [MP.dof(1:6).dof_r];
        A1(ldof,dof) = G.L;
        A1(dof,ldof) = transpose(G.L);
    end

end