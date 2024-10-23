function [A0,G,El,FLD,OHC,Ce,Ge,Ie,V0] = assemble_A0(Nd,El,OHC,IHC,FLD,MP)

    if MP.Visc == 0
        [A0,G,El,FLD,OHC,Ce,Ge,Ie,V0] = assemble_inviscid(Nd,El,OHC,IHC,FLD,MP);
    end

end

% #########################################################################
% ###################### INVISCID MATRICES ASSEMBLY #######################
% #########################################################################

function [A0,G,El,FLD,OHC,Ce,Ge,Ie,V0] = assemble_inviscid(Nd,El,OHC,IHC,FLD,MP)
% %
% % Assembe frequency-independent submatrices that do not need to update with stimulating frequencies
% %

    opt_cfld = MP.opt_cfld;

    %reduced degrees of freedom
    pdof = FLD.rdof;
    udof = Nd.rdof;
    edof = MP.dof(3).n;
    adof = MP.dof(4).n;
    odof = MP.dof(5).n;

    G.Apu = sparse(pdof,udof);
    G.Ape = sparse(pdof,edof);
    G.Apa = sparse(pdof,adof);
    G.Apo = sparse(pdof,odof);

    G.Auu = sparse(udof,udof);
    G.Aua = sparse(udof,adof);

    G.Aep = sparse(edof,pdof);
    G.Aeu = sparse(edof,udof);
    G.Aee = sparse(edof,edof);
    G.Aea = sparse(edof,adof);
    G.Aeo = sparse(edof,odof);

    G.Aap = sparse(adof,pdof);
    G.Aae = sparse(adof,edof);
    G.Aaa = sparse(adof,adof);
    G.Aao = sparse(adof,odof);

    G.Aop = sparse(odof,pdof);
    G.Aou = sparse(odof,udof);
    G.Aoe = sparse(odof,edof);
    G.Aoa = sparse(odof,adof);
    G.Aoo = sparse(odof,odof); 
    
    if opt_cfld == 1

        cdof = FLD.CFLD.rdof;

        G.Qpc = sparse(pdof,cdof);
        G.Auc = sparse(udof,cdof); % place holder, should be defined
        G.Acc = sparse(cdof,cdof); % place holder, should be defined

        G.Acu = transpose(G.Auc);        
        G.Aec = sparse(edof,cdof);
        G.Aac = sparse(adof,cdof);
        G.Aoc = sparse(odof,cdof);
    end

    a = attenuation_function(MP,Nd.Z(:),6,Nd.tdof,Nd.BC,1,3,10);
    El.mat(El.N) = struct('m',[],'k',[],'c',[]);
    [M,El] = createM(Nd,El);
    [K,El] = createK(Nd,El);
    [C,El] = createC2(Nd,El,MP);
    
    G.M = M;
    G.C = C;
    G.K = K;

    [G, FLD] = make_FE_FLDmatrix(MP,G,FLD,opt_cfld);
    
    [G,FLD] = make_AMatrix_FSI_2(G,Nd,El,FLD,MP);
    
    if opt_cfld == 1
        G = make_AMatrix_FSI_CFLD(Nd,FLD,G,MP);
        G = make_AMatrix_FLD_CFLD(MP,FLD,G); % permeability relation between the two fluid spaces
    end

    % C and G matrices at equilibrium state
    [Ce,Ge,Ie,V0] = mElectrical2_cmod_MSA(MP,OHC,IHC); 

    HB = [OHC.HB];
    G.Aau = make_Aau(MP,Nd,El,OHC);
    [G.Aaa,G.Iaa] = make_Aaa(MP,HB);       

    [G.Auo, ~] = make_Auo(MP,Nd,El,OHC,HB);
    [G.Aue, ~] = make_Aue(MP,Nd,El,OHC);
    % active comes later, because G.Auo and G.Aue are needed for rhs vector
    % for the zero frequency where po0 and Vm0 are non-zero

    G = make_Lagrange_constraints_matrices(G,MP,Nd,El);

    if opt_cfld == 1
        A0 = [ G.App, G.Apu, G.Ape, G.Apa, G.Apo, G.Qpc;
               G.Aup, G.Auu, G.Aue, G.Aua, G.Auo, G.Auc;    
               G.Aep, G.Aeu, G.Aee, G.Aea, G.Aeo, G.Aec;
               G.Aap, G.Aau, G.Aae, G.Aaa, G.Aao, G.Aac;
               G.Aop, G.Aou, G.Aoe, G.Aoa, G.Aoo, G.Aoc;
              G.Qpc', G.Acu,G.Aec',G.Aac',G.Aoc', G.Acc];
        ldof = length(MP.dof(7).dof_r);
        dof = length([MP.dof(1:6).dof_r]);
        LA = sparse(ldof,dof);
        LL = sparse(ldof,ldof);
        A0 = [ A0, LA.';
               LA, LL  ];
    else
        A0 = [ G.App, G.Apu, G.Ape, G.Apa, G.Apo;
               G.Aup, G.Auu, G.Aue, G.Aua, G.Auo;    
               G.Aep, G.Aeu, G.Aee, G.Aea, G.Aeo;
               G.Aap, G.Aau, G.Aae, G.Aaa, G.Aao;
               G.Aop, G.Aou, G.Aoe, G.Aoa, G.Aoo];
    end

    % structural Lagrange constraints
    ldof = MP.dof(8).n;    
    nd = Nd.Master_Slave(:,1);
    rdof = (nd-1)*6 + (1:3);
    rdof = Nd.mapF2R(rdof);
    row = logical(reshape(transpose(rdof),ldof,1));
    ldof_r = sum(row);
    dof = MP.rdof-ldof_r;
    col = false(1,dof);
    col(MP.dof(2).dof_r) = true;
    LA = sparse(ldof_r,dof);
    LA(:,col) = G.lu(row,:);
    LL = sparse(ldof_r,ldof_r);
    A0 = [ A0, LA.';
           LA, LL  ];

end

