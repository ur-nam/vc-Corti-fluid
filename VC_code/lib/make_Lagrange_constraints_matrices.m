function G = make_Lagrange_constraints_matrices(G,MP,Nd,El)

    NDOF = 6;
    udof = MP.dof(2).n;
    ldof = MP.dof(8).n;


    partsR = zeros(numel(Nd.Master_Slave)*3,4);
    cnt = 0;

    for ni = 1:length(Nd.Master_Slave)
        m = Nd.Master_Slave(ni,1);
        s = Nd.Master_Slave(ni,2);
        m_dof = (m-1)*NDOF + (1:3);
        s_dof = (s-1)*NDOF + (1:3);
        l_dof = (ni-1)*3 + (1:3);

        rdof = l_dof;
        cdof = [m_dof,s_dof];
        m = [1,0,0,-1,0,0;
             0,1,0,0,-1,0;
             0,0,1,0,0,-1;];
        [partsR,cnt] = fillR(rdof,cdof,m,partsR,cnt);
    end

    lu = sparse(partsR(:,1),partsR(:,2),partsR(:,3),ldof,udof);
    G.lu = lu(:,Nd.BC);

end