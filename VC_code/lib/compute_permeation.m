function vp = compute_permeation(MP, FLD, R, fi)

%     fprintf(1,'\n  :: Permeation for %4.2fkHz stim\n',R.freq(fi));
    CFLD = FLD.CFLD;

    kk = perm_scale(MP,CFLD.Nd.Z);
%     kk = kk.*CFLD.Nd.width;

    cdof = MP.dof(6).cdof;
    [~,idx] = sort(CFLD.Nd.x(CFLD.Nd.ind_Radi,1),'ascend');
    cpdof = ((1:CFLD.Nd.N)-1)*CFLD.Nd.NV + 3;
    cpdof_FSI = cpdof(CFLD.Nd.ind_Radi(idx));
    cc = R.Uf(cdof,fi);
    p_in = cc(cpdof_FSI);

    pdof = MP.dof(1).pdof;
    [~,idx] = sort(FLD.Nd.x(FLD.indBM,1),'ascend');
    spdof_FSI = FLD.indBM(idx);    
    pp = R.Uf(pdof,fi);
    p_out = pp(spdof_FSI);

    vp = kk.*(p_in - p_out);
end