function CFLD = set_CFLD_BC_1D(CFLD)

    Nd = CFLD.Nd;
    
    ind_Q = true(Nd.N,1);
%     ind_Q(1) = false;
    ind_Q(end) = false;

    ind_p = true(Nd.N,1);
    ind_p(1) = false;
%     ind_p(end) = false;
    CFLD.ind_Q = ind_Q;
    CFLD.ind_p = ind_p;
    CFLD.tdof = 2*Nd.N;
    CFLD.BC = [ind_Q.';ind_p.'];
    CFLD.BC = CFLD.BC(:);
    CFLD.rdof = sum(CFLD.BC);

end