function GV_nonlin = evaluate_GV_nonlin(MP,OHC,voltage,po,po0)

    nz = round(MP.length_BM/MP.dZ)+1;           % # of sections along longitudinal direction
    nv = 5;
    GV_nonlin = zeros(nv*nz,1);
    for iz = 1:nz
        rdof = (5*iz-4):5*iz;
        Gs = OHC(iz).count*OHC(iz).S.Gmax*(po(iz) - po0(iz));
        G = [Gs, -Gs, 0, 0, 0;
             -Gs, Gs, 0, 0, 0;
               0,  0, 0, 0, 0;
               0,  0, 0, 0, 0;
               0,  0, 0, 0, 0;];
        GV_nonlin(rdof,1) = GV_nonlin(rdof,1) + G*voltage(rdof,1);
    end

end