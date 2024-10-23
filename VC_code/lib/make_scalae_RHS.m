function b = make_scalae_RHS(MP,FLD,p_f)

    pdof = MP.dof(1).dof_r;

    b = zeros(MP.rdof,1);
    dtau_inv = 1/MP.dtau;
    
    row = pdof;
    b(row) = b(row) + dtau_inv*ones(FLD.rdof,1).*p_f(FLD.BC,1);
end