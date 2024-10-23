function d_stapes = compute_Stapes(MP,FLD,R)

    % non-zero frequency bins
    f_idx = find(R.freq ~= 0);
    freq = R.freq(f_idx);
    % approximated stapes displacement at each frequency
    y = 150; x = [0,10]; % reference location for pressure gradient
    [xx,yy] = meshgrid(x,y);
    idx = rangesearch(FLD.Nd.x,[5,150],50);
    d_stapes = zeros(length(f_idx),1);
    for fi = 1:length(f_idx)
        int_pre = scatteredInterpolant(FLD.Nd.x(idx{1},:), R.pre(idx{1},f_idx(fi)));
        pp = int_pre(xx,yy);
        dp = (pp(2)-pp(1))/(xx(2)-xx(1));
        ww = 2*pi*freq(fi);
        d_stapes(fi,1) = -dp(1)/((1j*ww)^2)/MP.rho*(1e3); % [nm]
    end

end