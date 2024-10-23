function post_MC(Nd,El,MP,fMP,FLD,OHC,R)

    loc = MP.loc + 1; % slit opening center

    idx = knnsearch(MP.xx',loc);

    % gain
    
    ref1 = 20*1e-6*10^(MP.Pstim/20); % [Pa]
    edof = MP.dof(3).edof;
    ee = R(2).Uf(edof,:);
    nv = 5;
    ref2 = 1e-3*ee((idx-1)*nv+1,:); % [V]

    yy1 = R(1).yBM(idx,:)./ref1; % [nm/Pa]
    yy2 = R(2).yBM(idx,:)./ref2; % [nm/V]

    figure(200);

    ax_tw_gain = subplot(2,1,1); % set(ax_tw,'Position',[0.12, 0.6, 0.3, 0.34]);
    hold(ax_tw_gain,'on');
    semilogy(ax_tw_gain,MP.freq,abs(yy1),'k','LineWidth',1);
    semilogy(ax_tw_gain,MP.freq,abs(yy2),'r','LineWidth',1);
    hold(ax_tw_gain,'off');

    ax_tw_phase = subplot(2,1,2);
    hold(ax_tw_phase,'on');
    semilogy(ax_tw_phase,MP.freq,unwrap(angle(yy1)),'k','LineWidth',1);
    semilogy(ax_tw_phase,MP.freq,unwrap(angle(yy2)),'r','LineWidth',1);
    hold(ax_tw_phase,'off');
end