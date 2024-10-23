function FLD = FLD_mesh(MP,H,FLD)

    H_05 = round(H/2);
    H_OCC = 25;
    L = MP.length_BM;
    h0 = 10;

    if MP.vMC == 1
        load('./VC_code/hinput/vMC_4mm.mat','p','t');
        Nd.x = p; El.node = t;
        El.type = 3;
        Lo = 1000;
        Ho = 80;
        Xo = 2500;
        Ls = 1000;
        Hs = 200;
        Xs = 2500;
        FLD.Lo = Lo;
        FLD.Xo = Xo;
        FLD.Ho = Ho;
        FLD.Ls = Ls;
        FLD.Xs = Xs;
        FLD.Hs = Hs;
    else
        % using quadriateral mesh
        load('./VC_code/hinput/mesh_12mm_cochlea_quad_4.mat','msh');
        Nd.x = msh.POS(:,[1,2]);    
        El.node = msh.QUADS(:,1:4);
        El.type = 4;     
    end

    El.N = size(El.node,1);
    Nd.N = size(Nd.x,1);
    FLD.Nd = Nd;
    FLD.El = El;

    FLD.H = H;
    FLD.h0 = h0;
    FLD.HOC = H_OCC;
    FLD.L = L;

end