function b = make_MET_RHS(opt_lin,G,wo,MP,V0,po0,Auo,Aue,po_f,voltage_f,po_f_rhs,GV_nonlin_f)

    udof = MP.dof(2).dof_r;
    edof = MP.dof(3).dof_r;
    odof = MP.dof(5).dof_r;

    b = zeros(MP.rdof,1); % reduced rhs vector

    dtau_inv = 1/MP.dtau;

    if ~opt_lin
        row = odof;
        b(row) = b(row) + po_f.*dtau_inv + po_f_rhs;
    
        row = edof;
        b(row) = b(row) + dtau_inv*G.Aee*voltage_f - GV_nonlin_f;
    end
    
    if wo == 0
        f_MET_rest = Auo*po0;
        row = udof;
        b(row) = b(row) + f_MET_rest;
        f_OHC_rest = Aue*V0;
        b(row) = b(row) + f_OHC_rest;
    end

end