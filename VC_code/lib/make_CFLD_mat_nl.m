function FLD = make_CFLD_mat_nl(fi, MP, FLD, Uf0, wk,opt_lin)
    
    MP.nonlin = ~opt_lin; % quick way to modify the code for new flag
    CFLD = FLD.CFLD;
    El = CFLD.El;
    Nd = CFLD.Nd;
    if MP.nonlin, convection = CFLD.convection; end
    Uf0 = reshape(Uf0,Nd.NV,Nd.N).';

    NV = Nd.NV; %number of variables
    N = El.type*NV; % number of dof for an element

    partsMA = zeros(N*N*El.N,4); cntA = 0;

    b = zeros(NV*Nd.N,1); % right hand side
    
    for ie = 1:El.N
        knodes = El.node(ie,:);
        if MP.nonlin
            [A_e,F_e] = element_matrix(ie,El,Nd,MP,convection(:,:,fi),Uf0(knodes,:), wk);
        else
            [A_e,F_e] = element_matrix(ie,El,Nd,MP,[],Uf0(knodes,:), wk);
        end
        [partsMA, cntA] = fillM(NV, knodes', A_e, partsMA, cntA);
        b = fillV(NV, knodes', F_e, b);
                 
    end

    Acc = sparse(partsMA(:,1), partsMA(:,2), partsMA(:,3) + 1j*partsMA(:,4), CFLD.tdof, CFLD.tdof);

    cdof = CFLD.BC;

%     zz = [fNd.Z.';fNd.Z.'];
%     a = attenuation_function(MP,zz(:),2,CFLD.tdof,CFLD.BC,1,3,10);
    CFLD.Acc = Acc(cdof,cdof);
    CFLD.b = b;
    FLD.CFLD = CFLD;
end

function [A_e,F_e] = element_matrix(ie,El,Nd,MP,convection,Uf, wk)

    [~, ~, w] = set_Gauss_local_variables(2,2,El.NQ);
    NQ = El.NQ;
    dtau_inv = 0;
    nu = MP.nu;
    NV = Nd.NV;
    rho_inv = 1/MP.rho;
    A_e = zeros(NV*El.type);
    F_e = zeros(NV*El.type,1);
    L = zeros(4,NV*El.type); %first dim is number of equations
    hs = El.elm_hs(:,:,ie);
    CONV = zeros(1,2*NQ);
    if MP.nonlin, CONV(1,:) = convection(ie,:,1); end
    for iq = 1:NQ
        
        cf = hs(iq)*w(iq);

        cpsi = El.elm_psi(:,iq,ie);
        cgpsix = El.elm_gpsi(1:El.type,iq,ie);
        cgpsiy = El.elm_gpsi(1 + El.type:end,iq,ie);

        u = cpsi'*Uf(:,1);
        v = cpsi'*Uf(:,2);        
        dpdx = cgpsix'*Uf(:,3);
        dpdy = cgpsiy'*Uf(:,3);
        dwdx = cgpsix'*Uf(:,4);
        dwdy = cgpsiy'*Uf(:,4);

        theta = 1;
        a_u = theta; a_v = theta; a_p = theta; a_w = theta; 
        A1 = [0 0 a_p*rho_inv 0;0 0 0 -a_w*nu; 0 -1 0 0; 1 0 0 0];
        A2 = [0 0 0 a_w*nu;0 0 a_p*rho_inv 0; 1 0 0 0; 0 1 0 0];
        A0 = [dtau_inv + a_u*1j*wk 0 0 0; 0 dtau_inv + a_v*1j*wk 0 0;...
            0 0 0 1; 0 0 0 0];
        f = [dtau_inv*u - ((1-a_u)*1j*wk*u + CONV(1,iq) + (1-a_p)*rho_inv*dpdx + (1-a_w)*nu*dwdy);...
             dtau_inv*v - ((1-a_v)*1j*wk*v + CONV(1,iq+NQ) + (1-a_p)*rho_inv*dpdy - (1-a_w)*nu*dwdx);...
             0; 0];
       
        for in = 1:El.type

            psi = El.elm_psi(in,iq,ie);
            gpsix = El.elm_gpsi(in,iq,ie);
            gpsiy = El.elm_gpsi(in + El.type,iq,ie);

            L(:,NV*(in-1)+1:NV*in) = gpsix*A1 + gpsiy*A2 + psi*A0;
            
        end
        
        A_e = A_e + (L.')*L*cf;
        F_e = F_e + (L.')*f*cf;
        
    end
end
