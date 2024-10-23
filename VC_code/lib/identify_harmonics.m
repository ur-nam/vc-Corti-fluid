function MP = identify_harmonics(MP,U,threshold)
    
    N = ndims(U); % this is to be able to take matrices that are of different dimensions
    yy = reshape(mean(abs(U),1:N-1),1,MP.L);

    % determine frequency indices
    signal_threshold = threshold; % dB
    ref = max(yy);
    
    ind = find((20*log10(yy/ref) > - signal_threshold) & MP.ffsym > 0);
        
    ind = [MP.sind,ind];
    ind = unique(ind);

    MP.ind = ind;
%     if length(MP.ind) > 40
%         MP.ind = MP.ind(1:40);
%     end
    figure(21); clf; ax = axes;
    plot_spectrum(ax,MP,yy);
end