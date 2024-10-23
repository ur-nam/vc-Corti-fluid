function [u,v,freq] = scala_fluid_motion(MP,FLD,R,nf)

    Nf = length(nf);

    % (N')*B*pre = -rho*(N')*N*a: weak form of dpdx = -rho*a
    [M, M_l, Dx, Dy, B] = make_matrices(FLD);

    freq = zeros(1,Nf);
    u = zeros(FLD.Nd.N,Nf);
    v = zeros(FLD.Nd.N,Nf);

    icount = 0;
    for fi = nf

        icount = icount + 1;

        freq(icount) = MP.freq(fi);

        pre = R.pre(:,fi);

%         [A,bx,by] = lsq_matrices(MP,FLD,pre,fi);
%         u(:,icount) = A\bx;
%         v(:,icount) = A\by;

        ww = 2*pi*MP.freq(fi);
        pre = R.pre(:,fi);
        b = -Dx*pre/FLD.rho;
        u(:,icount) = (1i*ww*M_l)\b;
%         u(:,icount) = b/(1i*ww);
        b = -Dy*pre/FLD.rho;
        v(:,icount) = (1i*ww*M_l)\b;
%         v(:,icount) = b/(1i*ww);

        % kinetic energy 

    end



end

function [A,bx,by] = lsq_matrices(MP,FLD,pre,fi)

    Nd = FLD.Nd;
    El = FLD.El;
    
    NV = 1;
    N = El.type*NV;
    NQ = El.NQ;

    ww = 2*pi*MP.freq(fi);

    partsA = zeros(N*N*El.N,4); cntA = 0;
    bx = zeros(NV*Nd.N,1); % right hand side
    by = zeros(NV*Nd.N,1); % right hand side

    [~, ~, w] = set_Gauss_local_variables(2,2,El.NQ);

    for ie = 1:El.N        

        knodes = El.node(ie,:);

        A_e = zeros(NV*El.type);
        Fx_e = zeros(NV*El.type,1);
        Fy_e = zeros(NV*El.type,1);
        L = zeros(1,NV*El.type); %first dim is number of equations
        
        hs = El.elm_hs(:,:,ie);

        for iq = 1:NQ
    
            cf = hs(iq)*w(iq);

            cgpsix = El.elm_gpsi(1:El.type,iq,ie);
            cgpsiy = El.elm_gpsi(1 + El.type:end,iq,ie);    
    
            dpdx = cgpsix'*pre(knodes,1);
            dpdy = cgpsiy'*pre(knodes,1);

            A0 = -FLD.rho*1i*ww;
            fx = dpdx;
            fy = dpdy;

            for in = 1:El.type
                psi = El.elm_psi(in,iq,ie);
                L(:,NV*(in-1)+1:NV*in) = psi*A0;
            end

            A_e = A_e + (L.')*L*cf;
            Fx_e = Fx_e + (L.')*fx*cf;
            Fy_e = Fy_e + (L.')*fy*cf;
        end

        [partsA, cntA] = fillM(NV, knodes', A_e, partsA, cntA);
        bx = fillV(NV, knodes', Fx_e, bx);
        by = fillV(NV, knodes', Fy_e, by);

    end

    A = sparse(partsA(:,1), partsA(:,2), partsA(:,3) + 1j*partsA(:,4), FLD.tdof, FLD.tdof);

end

function [M, M_l, Dx, Dy, B] = make_matrices(FLD)

    Nd = FLD.Nd;
    El = FLD.El;
    
    NV = 1;
    N = El.type*NV;

    partsM = zeros(N*N*El.N,4); cntM = 0;
    partsM_l = zeros(N*N*El.N,4); cntM_l = 0;
    partsDx = zeros(N*N*El.N,4); cntDx = 0;
    partsDy = zeros(N*N*El.N,4); cntDy = 0;
    partsB = zeros(N*N*El.N,4); cntB = 0;

    for ie = 1:El.N
        knodes = El.node(ie,:);
        [Me, Me_l, De_x,De_y, Be] = element_matrix(ie,Nd,El,NV);
        [partsM, cntM] = fillM(NV, knodes', Me, partsM, cntM);
        [partsM_l, cntM_l] = fillM(NV, knodes', Me_l, partsM_l, cntM_l);
        [partsDx, cntDx] = fillM(NV, knodes', De_x, partsDx, cntDx);
        [partsDy, cntDy] = fillM(NV, knodes', De_y, partsDy, cntDy);
        [partsB, cntB] = fillM(NV, knodes', Be, partsB, cntB);
    end

    M = sparse(partsM(:,1), partsM(:,2), partsM(:,3) + 1j*partsM(:,4), FLD.tdof, FLD.tdof);
    M_l = sparse(partsM_l(:,1), partsM_l(:,2), partsM_l(:,3) + 1j*partsM_l(:,4), FLD.tdof, FLD.tdof);
    Dx = sparse(partsDx(:,1), partsDx(:,2), partsDx(:,3) + 1j*partsDx(:,4), FLD.tdof, FLD.tdof);
    Dy = sparse(partsDy(:,1), partsDy(:,2), partsDy(:,3) + 1j*partsDy(:,4), FLD.tdof, FLD.tdof);
    B = sparse(partsB(:,1), partsB(:,2), partsB(:,3) + 1j*partsB(:,4), FLD.tdof, FLD.tdof);

end

function [M, M_l, Dx, Dy, B] = element_matrix(ie,Nd,El,NV)

    NQ = El.NQ;
    [~, ~, w] = set_Gauss_local_variables(2,2,NQ);

    M = zeros(NV*El.type);
    Dx = zeros(NV*El.type);
    Dy = zeros(NV*El.type);
    B = zeros(NV*El.type);
    hs = El.elm_hs(:,:,ie);

    for iq = 1:NQ

        cf = hs(iq)*w(iq);

        cpsi = El.elm_psi(:,iq,ie); % [N,1]
        cgpsix = El.elm_gpsi(1:El.type,iq,ie); % [N,1]
        cgpsiy = El.elm_gpsi(1 + El.type:end,iq,ie); % [N,1]

        M = M + cf*cpsi*(cpsi.');
        Dx = Dx + cf*cpsi*(cgpsix.');
        Dy = Dy + cf*cpsi*(cgpsiy.');
        B = B + cf*(cgpsix*(cgpsix.') + cgpsiy*(cgpsiy.')); % laplacian

    end

    M_l = diag(sum(M,1));

end