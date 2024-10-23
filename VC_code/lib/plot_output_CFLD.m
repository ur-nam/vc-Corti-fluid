function plot_output_CFLD(Uf,fMP,MP,Nd,FLD)

    wk = 2*pi*fMP.freq;
    xx = MP.xx;
    % pressure
    cdof = reshape(MP.dof(6).cdof,2,[]);
    ind = find(fMP.ffsym == fMP.freq);
    cp = Uf(cdof(2,:),ind);

    bcTM = FLD.indTM;
    [~,idx] = sort(FLD.Nd.x(bcTM,1),'ascend');
    bcTM = bcTM(idx);
    bcBM = FLD.indBM;
    [~,idx] = sort(FLD.Nd.x(bcBM,1),'ascend');
    bcBM = bcBM(idx);

    pdof = MP.dof(1).pdof;
    pre = Uf(pdof,ind);
    tp = pre(bcTM);
    bp = pre(bcBM);
    cmap = lines(3);
    figure(44); clf;
    subplot(3,1,1);
    plot(xx,real(cp),'Color',cmap(1,:),'LineStyle','-'); hold on
    plot(xx,real(tp),'Color',cmap(2,:),'LineStyle','-');
    plot(xx,real(bp),'Color',cmap(3,:),'LineStyle','-');
    plot(xx,abs(cp),'Color',cmap(1,:),'LineStyle','--');
    plot(xx,abs(tp),'Color',cmap(2,:),'LineStyle','--');
    plot(xx,abs(bp),'Color',cmap(3,:),'LineStyle','--');
    legend({'cp [Pa]','tp [Pa]','bp [Pa]'});
    xlim([min(MP.xx),max(MP.xx)]);
    hold off
    % area change and flux
    
    CFLD = FLD.CFLD;
    dA = zeros(CFLD.Nd.N,1);

    ndof = 6;
    udof = MP.dof(2).udof;
    uu_f = Uf(udof,ind);
    idx = reshape(transpose((1:Nd.tdof)),ndof,Nd.N);
    dx_f = uu_f(idx(1,:));
    dy_f = uu_f(idx(2,:));

    for in = 1:CFLD.Nd.N
        knodes = CFLD.Nd.perim(in,:);
        dx = dx_f(knodes);
        dy = dy_f(knodes);
        H = CFLD.Nd.H(in,:);
        dA(in) = 1j*wk*dot(H,[dx;dy]);
    end

    cq = Uf(cdof(1,:),ind);
    dcqdz = gradient(cq,xx*1e3);

    kk = MP.kk;
    width = FLD.CFLD.Nd.width;
    dP = cp - bp;
    chi = MP.dZ*kk*width(:).*dP(:);

    cmap = lines(3);
    subplot(3,1,2);
    plot(xx,real(dcqdz),'Color',cmap(1,:),'LineStyle','-'); hold on
    plot(xx,real(dA),'Color',cmap(2,:),'LineStyle','-');
    plot(xx,real(chi),'Color',cmap(3,:),'LineStyle','-');
    plot(xx,abs(dcqdz),'Color',cmap(1,:),'LineStyle','--');
    plot(xx,abs(dA),'Color',cmap(2,:),'LineStyle','--');
    plot(xx,abs(chi),'Color',cmap(3,:),'LineStyle','--');
    legend({'dcqdz [\mum^2/ms]','dAdt [\mum^2/ms]','permeation [\mum^2/ms]'});
    xlim([min(MP.xx),max(MP.xx)]);
    hold off

    % permeation

    kk = MP.kk;
    dy = dy_f(FLD.nBot);
    dP = cp - bp;
    vs = 1j*wk*dy;
    vp = -kk*width(:).*dP(:);

    cmap = lines(2);
    subplot(3,1,3);
    plot(xx,real(vs),'Color',cmap(1,:),'LineStyle','-'); hold on
    plot(xx,real(vp),'Color',cmap(2,:),'LineStyle','-');
    plot(xx,abs(vs),'Color',cmap(1,:),'LineStyle','--');
    plot(xx,abs(vp),'Color',cmap(2,:),'LineStyle','--');
    legend({'vs [\mum/ms]','vp [\mum/ms]'});
    xlim([min(MP.xx),max(MP.xx)]);
    hold off
    drawnow
end