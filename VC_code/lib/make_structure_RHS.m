function b = make_structure_RHS(MP,Nd,uu_f)

    udof = MP.dof(2).dof_r;
    
    b = zeros(MP.rdof,1); % reduced rhs vector
    dtau_inv = 1/MP.dtau;

    row = udof;
    b(row) = b(row) + dtau_inv*ones(Nd.rdof,1).*uu_f(Nd.BC,1);

end