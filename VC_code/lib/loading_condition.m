function [F, Uf, FLD] = loading_condition(MP,fMP,FLD,V_lin,Ie,po0,Uf)

    F = zeros(MP.tdof,fMP.L,'like',1j);

    pre_in = zeros(FLD.tdof,fMP.L);
    
    if MP.vMC == 1
        ind_pre = strncmp(FLD.name,'BCS',3);
    else
        ind_pre = strncmp(FLD.name,'OW',2);
    end

    col = knnsearch(fMP.ffsym.',fMP.freq).';
    
    if length(MP.Pstim) > 1
        amplitude = MP.Pstim/20;
    else
        stim_f = zeros(1,fMP.L);
        col = knnsearch(fMP.ffsym.',fMP.freq).';
        amplitude = ones(1,fMP.bins-1); % minus the DC component
        stim_f(col) = amplitude;
        stim_f(:,fMP.L-col+2) = conj(amplitude);
        stim = inverse_fourier_transform(stim_f);
        factor = 1/sqrt(mean(stim.^2));
        amplitude = factor*amplitude*MP.Pstim/20;
    end

    if (MP.Estim == 1)
        pre_in(ind_pre,col) = 0;
    else
        pre_in(ind_pre,col) = ones(sum(ind_pre),1)*20*1e-6*10.^amplitude*exp(1j*eps); % prescribed pressure applied
    end

    pre_in(ind_pre,fMP.L-col+2) = conj(pre_in(ind_pre,col));
    RHS = -FLD.L*pre_in;

    row = MP.dof(1).pdof;
    F(row,:) = F(row,:) + RHS;

    if MP.key_fOHC || MP.key_fMET % active cochlea condition
        col = (fMP.ffsym == 0);
        row = MP.dof(3).edof;
        F(row,col) = Ie + V_lin;
    end
    
    if MP.Estim == 1
        nv = 5; % 5 circuit nodes
        I_in = zeros(MP.dof(3).n,fMP.L);
        col = knnsearch(fMP.ffsym.',fMP.freq).';
        amp = MP.Istim*1e9; % [mA] to [pA]
        ind_StV = 1:nv:MP.dof(3).n; % current nodes for StV
        I_in(ind_StV,col) = amp;
        I_in(ind_StV,fMP.L-col+2) = conj(I_in(ind_StV,col));
        row = MP.dof(3).edof;
        F(row,:) = F(row,:) + I_in;
    end

    if ~MP.nonlin
        col = (fMP.ffsym == 0);
        row = MP.dof(5).odof;
        F(row,col) = po0;
    end

    F = F(MP.BC,:);

    FLD.pre_in = pre_in;

    row = MP.dof(1).pdof;
    Uf(row(ind_pre),:) = pre_in(ind_pre,:);
    
end