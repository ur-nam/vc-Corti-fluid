function L = make_peristalsis_effect(MP, FLD, wk)

    ndof = 6;
    CFLD = FLD.CFLD;
    Nd = CFLD.Nd;
    [~,idx] = sort(Nd.x(Nd.ind_Radi,1),'ascend');
    vdof = ((1:Nd.N)-1)*Nd.NV + 2; %transverse velocity dof of Corti
    vdof_FSI = vdof(Nd.ind_Radi(idx));
    cpdof = ((1:Nd.N)-1)*Nd.NV + 3;
    cpdof_FSI = cpdof(Nd.ind_Radi(idx));
    [~,idx] = sort(FLD.Nd.x(FLD.indBM,1),'ascend');
    spdof_FSI = FLD.indBM(idx);   

    kk = perm_scale(MP,Nd.Z);
    kk = kk.*Nd.width;

    % dA/dt -2*pi*r*(1j*w*r) = 0;
    % dA/dt - 2*pi*r*vs = 0;
    % vf = vs + vp;
    % dA/dt - 2*pi*r*vf + psi = 0;
    % dA/dt = 1j*w*H;
    % psi = kk(CP - SP);
    % dA/dt - 2*pi*r*vf + kk*CP - kk*SP = 0;
    p_num = Nd.p_num;
    d_num = 2*p_num;

    partsL = zeros(Nd.nz*(d_num+3),4); cntL = 0;
    for i = 1:Nd.nz
        r = Nd.radi(i);
        k = kk(i);% *exp(-MP.p_beta*wk);
        H = Nd.H(i,:);
        H = [H(1:p_num);H(p_num+1:2*p_num)]; H = H(:); % [dx1,dy1,dx2,dy2,...]
        m = [1j*wk*H;-2*pi*r;-k;k];

        perim = Nd.perim(i,:);
        sdof = [1;2] + (perim - 1)*ndof; % structure dof
        col = [MP.dof(2).udof(sdof(:)),MP.dof(6).cdof(vdof_FSI(i)),MP.dof(6).cdof(cpdof_FSI(i)),MP.dof(1).pdof(spdof_FSI(i))];
        row = i*ones(size(col));
        [partsL,cntL] = fillT(row,col,m,partsL,cntL);
    end
    idx = 1:sum([MP.dof(1:6).n]);
    BC_s = MP.BC(idx,1);
    L = sparse(partsL(:,1), partsL(:,2), partsL(:,3) + 1j*partsL(:,4), CFLD.Nd.nz, sum([MP.dof(1:6).n]));
    BC_l = reshape(FLD.CFLD.BC,4,FLD.CFLD.Nd.N);
    BC_l = BC_l(2,FLD.CFLD.Nd.ind_Radi);
    L = L(BC_l,BC_s);
end