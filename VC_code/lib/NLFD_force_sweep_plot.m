function NLFD_force_sweep_plot

    n = 601; % 6mm half length. can be changed
    rdir = './results/force_sweep_NLFD/';
    list = ls(rdir);

    icount = 0;
    for i = 1:size(list,1)
        entry = list(i,:);
        if ~contains(entry, 'NLFD')
            continue;
        end
        icount = icount + 1;
        load([rdir,entry]);
        f(icount) = MP.stim;
        uf = Uf(MP.udof,:);        
        % reconstruct wave using frequency information
        xfHB = uf(n,:);
        tt = linspace(0,MP.tt(end),201);
        yy = real(xfHB(1,MP.ind)*exp(1i*2*pi*(MP.ffsym(MP.ind).').*tt));
        xHB(icount) = max(yy);
        % reconstruct open probability wave using frequency information
        ef = Uf(MP.edof,:);
        pf = ef(n,:);
        yy = real(pf(1,MP.ind)*exp(1i*2*pi*(MP.ffsym(MP.ind).').*tt));
        po(icount) = max(yy);
    end

    figure;
    subplot(1,2,1);
    plot(xHB*1e3,f);
    subplot(1,2,2);
    plot(xHB*1e3,po);

end