function plot_results(ax,MP,FLD,R)

    % traveling wave peaks
    [~,peak_idx] = max(abs(R.yBM),[],1);

    % stapes displacement at each frequency
    y = 150; x = [0,10];
    [xx,yy] = meshgrid(x,y);
    idx = rangesearch(FLD.Nd.x,[5,150],50);
    d_stapes = zeros(MP.Nf,1);
    for fi = 1:MP.Nf
        int_pre = scatteredInterpolant(FLD.Nd.x(idx{1},:), R.pre(idx{1},fi));
        pp = int_pre(xx,yy);
        dp = (pp(2)-pp(1))/(xx(2)-xx(1));
        ww = 2*pi*R.freq(fi);
        d_stapes(fi,1) = -dp(1)/((1j*ww)^2)/MP.rho*(1e3); % [nm]
    end

    % frequency location map

    axes(ax(1)); cla;
    xx = MP.xx(peak_idx);
    semilogy(xx,R.freq,'-k','Marker','.','LineWidth',1.0), hold on;
    semilogy(MP.xx,loc2freq(MP.xx),'-r','linewidth',1.0);
    title('Frequency-location rel.');
    xlabel('X (mm)','Fontname','Arial','Fontsize',12);
    ylabel('Frequency (kHz)','Fontname','Arial','Fontsize',12);
    set(gca,'xtick',0:2:12,'xticklabel',{'0','','4','','8','','12'},'xlim',[min(MP.xx),max(MP.xx)]);
    set(gca,'ylim',[0.2,60],'ytick',[0.3, 1, 3, 10, 30]); ylim([0.3 30]); xlim([0 12]);
    grid on;
    naxis;

    % traveling wave
    axes(ax(2)); cla;
    loc = [2,6,10];
    nn = length(loc);
    poi_idx = knnsearch(xx.',loc.');
    cmap = lines(3);
    for i = 1:nn
        b = knnsearch(MP.xx(:),freq2loc(loc2freq(xx(poi_idx(i)))*5.6569));
        a = knnsearch(MP.xx(:),freq2loc(loc2freq(xx(poi_idx(i)))*0.3536));
        yy1 = abs(R.yBM(:,poi_idx(i))/d_stapes(poi_idx(i)));
        yy2 = abs(R.yTM(:,poi_idx(i))/d_stapes(poi_idx(i)));
%         yy3 = abs(R.aDC(:,f_idx(poi_idx(i)))/d_stapes(poi_idx(i)));
        semilogy(MP.xx(b:a),yy1(b:a),'Color',cmap(1,:),'LineWidth',1.0); hold on
        semilogy(MP.xx(b:a),yy2(b:a),'Color',cmap(2,:),'LineWidth',1.0);
%         semilogy(MP.xx,yy3,'Color',cmap(3,:),'LineWidth',1.0);
    end
    xlabel('X (mm)','Fontname','Arial','Fontsize',12);
    set(gca,'xtick',0:2:12,'xticklabel',{'0','','4','','8','','12'},'xlim',[min(MP.xx),max(MP.xx)]);
    grid on;
    ylim([1e-1 1e3])
    xlim([0 12])
    naxis
end