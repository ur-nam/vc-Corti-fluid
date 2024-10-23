function plot_spectrum(ax,fMP,uu)    
    yy = mean(abs(uu),1);
    semilogy(ax,fMP.ffsym,yy,'k-');
    hold on;
    semilogy(ax,fMP.ffsym(fMP.ind),yy(:,fMP.ind),'rx');
    ylim([eps+mean(yy)*0.3, max(yy)*3])
    drawnow;
end