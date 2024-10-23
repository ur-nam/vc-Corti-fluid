function plot_cross_section(MP,Nd,El,R,loc,nf)
% loc: location
% nf: # reference of the frequency

    figure; clf; ax = axes;
    NDOF = 6;
    pfact = 30;

    nz = knnsearch(MP.xx(:),loc);
    nn = ((1:Nd.Nc)-1)*Nd.Nr + nz;

    xx = Nd.X(nz,:);
    yy = Nd.Y(nz,:);

    udof = MP.dof(2).udof;
    uu = reshape(R.Uf(udof,nf),NDOF,Nd.N);
    dx = uu(1,nn);
    dy = uu(2,nn);
    dz = uu(3,nn);
    
    ps = pfact/max([abs(uu(1,:)),abs(uu(2,:))]);

    hold(ax,'on');
    n_sec = find(abs(Nd.Z + MP.loc*1e3 - loc*1e3) < 5);
    e_idx = find(ismember(El.Nd1,n_sec));

    for ie=1:El.N   
        if ismember(ie,e_idx)
        nds = [El.Nd1(ie), El.Nd2(ie)];               
            %undeformed shadow
            plot3(ax,Nd.X(nds), Nd.Y(nds), Nd.Z(nds),'Color',[0.8 0.8 0.8],'LineWidth', 1.0);
            [R,T] = rotation_matrix(El,ie);
            d = uu(:,nds);
            coord = [Nd.X(nds);Nd.Y(nds);Nd.Z(nds)];
            d = R*d(:); % move displacement values to local frame
            L = El.L(ie);
            nl = 20;
            x = linspace(0,El.L(ie),nl);
            xx_l = zeros(3,nl);
            dy_l = zeros(1,nl);
            dx_l = zeros(1,nl);
            dz_l = zeros(1,nl);
            O = zeros(1,nl);
            for ii = 1:numel(x)
                % this is in x-y plane
                N = [(L-x(ii))/L, x(ii)/L]; %[u1,u2]
                xx_l(:,ii) = coord*N.';
                dx_l(ii) = real(N*d([1;7]));
                N = [1-3*x(ii)^2/L^2+2*x(ii)^3/L^3, x(ii)-2*x(ii)^2/L+x(ii)^3/L^2,...
                     3*x(ii)^2/L^2-2*x(ii)^3/L^3,-x(ii)^2/L+x(ii)^3/L^2]; %[v1,t1,v2,t2] local
                dy_l(ii) = real(N*d([2;6;8;12]));
                dz_l(ii) = real(N*d([3;5;9;11]));
            end
            dd = T.'*ps*[dx_l;dy_l;dz_l];
            plot3(ax,xx_l(1,:) + dd(1,:), xx_l(2,:) + dd(2,:), xx_l(3,:) + dd(3,:),'Color',[0,0,0],'LineWidth', 2.5);
        end
    end

%     plot3(ax,xx + ps*real(dx), yy + ps*real(dy), 'o','Color','none','MarkerFaceCol',[0, 0.4470, 0.7410],'MarkerSize',4);

%     tdof = length(Nd.mapF2R);
%     rdof = length(F);
%     tmp = zeros(tdof,1);
%     tmp(Nd.mapR2F(1:rdof)) = F;
%     F = 3e-2*reshape(tmp,6,Nd.N);
%     fx = F(1,nn);
%     fy = F(2,nn);
% 
%     quiver(ax,xx + ps*dx, yy + ps*dy, fx, fy,'LineWidth',2,...,
%         'Autoscale','off','Color','r');
%     daspect(ax,[1 1 1]);
%     maxF = max(abs(F));

    padding = 0.1*(max(xx) - min(xx));
    xlim([min(xx) - padding, max(xx) + padding]); % should be adjusted manually for now
    ylim([min(yy) - padding, max(yy) + padding]); % should be adjusted manually for now

    view(135,-45)
    daspect(ax,[1,1,1])
end

function  [R,T] = rotation_matrix(El,ie)
    
    if abs(El.dir(ie,3)-1)<0.1,
        ez = [1, 0, 0];
        % note that most elements alligned either radial(x) direction
        % or longitudianl(z) direction
    else
        ez = [0, 0, 1];
    end

    xv = El.dir(ie,:);

    yv = mycross(ez,xv);
    yv = yv/norm(yv);
    zv = mycross(xv,yv);
    zv = zv/norm(zv);

    T = [xv;yv;zv];
% %     R = zeros(9);
% %     for ii=1:4,
% %         R((ii-1)*3+1:(ii-1)*3+3,(ii-1)*3+1:(ii-1)*3+3) = T;
% %     end
    O33 = zeros(3);    
    R = [T, O33, O33, O33;  O33, T, O33, O33;  O33, O33, T, O33;   O33, O33, O33, T;];
end