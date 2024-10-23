function h = initiate_plots

    h = gobjects(7,1);
    cmap = lines(2);
    labels = {'yBM','po'};
    figure(31);
    ax = gca;
    h(1) = animatedline(ax,Color=cmap(1,:),Marker='.',MaximumNumPoints=100);
    set(gca,'XScale','log','YScale','log'); set(gca, 'YColor', cmap(1,:));
    yyaxis right;
    h(2) = animatedline(ax,Color=cmap(2,:),Marker='.',MaximumNumPoints=100);
    set(gca,'XScale','log','YScale','log'); set(gca, 'YColor', cmap(2,:));
    legend(ax,labels);

%     % two-tone suppression
%     cmap = lines(3);
%     labels = {'(f2-f1)/f2','(f2+f1)f2','(2f2-f1)/f2'};
%     figure(32);
%     ax = gca;
%     h(3) = animatedline(ax,Color=cmap(1,:),Marker='.',MaximumNumPoints=100);
%     set(gca,'YScale','linear');
%     h(4) = animatedline(ax,Color=cmap(2,:),Marker='.',MaximumNumPoints=100);
%     set(gca,'YScale','linear');
%     h(5) = animatedline(ax,Color=cmap(3,:),Marker='.',MaximumNumPoints=100);
%     set(gca,'YScale','linear');
%     legend(ax,labels);

    h(6) = figure(33);
    h(7) = figure(34); 

end