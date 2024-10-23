function dA = compute_area_change(MP, Nd, FLD, R, fi)

%     fprintf(1,'\n  :: Area change for %4.2fkHz stim\n',R.freq(fi));
    ndof = 6;
    CFLD = FLD.CFLD;
    udof = MP.dof(2).udof;
    uu = reshape(R.Uf(udof,fi),ndof,Nd.N);    
    p_num = CFLD.Nd.p_num;

    dA = zeros(1,CFLD.Nd.nz);
    for i = 1:CFLD.Nd.nz
        H = CFLD.Nd.H(i,:);
        H = [H(1:p_num);H(p_num+1:2*p_num)]; H = H(:); % [dx1,dy1,dx2,dy2,...]
        d = [uu(1,CFLD.Nd.perim(i,:));uu(2,CFLD.Nd.perim(i,:))];
        dA(i) = transpose(H)*d(:);
    end

end

