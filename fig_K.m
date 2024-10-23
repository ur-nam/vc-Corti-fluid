function fig_K(r_opt)

    % r_opt: restart option; 1:restart, 0:use saved results
    if ~exist('r_opt','var')
        r_opt = 0;
    end

    add_paths
    if r_opt
        z={'./VC_code/hinput/','./VC_code/houtput/','6','12','[0.5,25]',...
            '11','20','[1e3,1e0]','0','1','1','0','0.1'};
        load('./VC_code/hinput/m_coef.mat','coef');
        set_Global_coefficients(coef);
        [Nd,El,MP,fMP,FLD,OHC,R] = NLFD_virtualCochleaHarmonic(z);
    else
        load('./VC_code/houtput/081023_2332_020.mat','MP','Nd','El','R');
    end
    plot_fig(MP,Nd,El,R);
    
end

function plot_fig(MP,Nd,El,R)

    freq = MP.freq;
    xx = MP.xx;

    % Nakajima 2009 intracochlear pressure
    p_ref = (20e-6)*(10^((MP.Pstim)/20));
    Z_c = 20e9; % Gohm
    footplate_a = 3.2e-6; % m2
    d_ref = 1e3*(p_ref/Z_c)/footplate_a./(2*pi*MP.freq); % stapes displacement

    NDOF = 6;
    udof = MP.dof(2).udof;
    udof = reshape(udof,NDOF,Nd.N);
    
    nBM = strcmp(Nd.name,'AP');
    nRL = strcmp(Nd.name,'BB');
    nDC = strcmp(Nd.name,'DD');

    figure; clf; set(gcf,'Position',[300, 60, 800, 900]);
    h = gobjects(13,1);
    
    pfact = 0.25;
    % Active

    xBM = R(2).Uf(udof(1,nBM),:);
    yBM = R(2).Uf(udof(2,nBM),:);
    zBM = R(2).Uf(udof(3,nBM),:);
    xRL = R(2).Uf(udof(1,nRL),:);
    yRL = R(2).Uf(udof(2,nRL),:);
    zRL = R(2).Uf(udof(3,nRL),:);
    xDC = R(2).Uf(udof(1,nDC),:);
    yDC = R(2).Uf(udof(2,nDC),:);
    zDC = R(2).Uf(udof(3,nDC),:);

    xBMp = R(1).Uf(udof(1,nBM),:);
    yBMp = R(1).Uf(udof(2,nBM),:);
    zBMp = R(1).Uf(udof(3,nBM),:);
    xRLp = R(1).Uf(udof(1,nRL),:);
    yRLp = R(1).Uf(udof(2,nRL),:);
    zRLp = R(1).Uf(udof(3,nRL),:);
    xDCp = R(1).Uf(udof(1,nDC),:);
    yDCp = R(1).Uf(udof(2,nDC),:);
    zDCp = R(1).Uf(udof(3,nDC),:);
    
    loc = 401;
    [~,peak_idx] = max(abs(yBM(loc,:)));
    
    f_tail = knnsearch(freq',freq(peak_idx)*0.4);
    f_peak = knnsearch(freq',freq(peak_idx)*0.9);

    b = knnsearch(freq',freq(peak_idx)*1.5);
    a = knnsearch(freq',freq(peak_idx)*0.25);
    
    n = 1;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.1125,0.790,0.2210,0.1600]);
    xlim(h(n),[-2.25 1])
    xticks(h(n),-2:0.5:1)
    xticklabels(h(n),{'','','','','','',''})
    ylim(h(n),[-10 60])
    yticks(h(n),0:20:60)
    yticklabels(h(n),{'0','20','40','60'})
    ylabel(h(n),'Gain re. stapes (dB)')
    title(h(n),'Transverse','FontWeight','normal','FontSize',11);

    n = 2;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.3910,0.790,0.2210,0.1600]);
    xlim(h(n),[-2.25 1])
    xticks(h(n),-2:0.5:1)
    xticklabels(h(n),{'','','','','','',''})
    ylim(h(n),[-10 60])
    yticks(h(n),0:20:60)
    yticklabels(h(n),{'','','',''})
    title(h(n),'Radial','FontWeight','normal','FontSize',11);

    n = 3;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.6660,0.790,0.2210,0.1600]);
    xlim(h(n),[-2.25 1])
    xticks(h(n),-2:0.5:1)
    xticklabels(h(n),{'','','','','','',''})
    ylim(h(n),[-10 60])
    yticks(h(n),0:20:60)
    yticklabels(h(n),{'','','',''})
    title(h(n),'Longitudinal','FontWeight','normal','FontSize',11);

    n = 4;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.1125,0.5930,0.2210,0.1600]);
    xlim(h(n),[-2.25 1])
    xticks(h(n),-2:0.5:1)
    xticklabels(h(n),{'','','','','','',''})
    ylim(h(n),[-20 40])
    yticks(h(n),-20:20:40)
    yticklabels(h(n),{'-20','0','20','40'})
    ylabel(h(n),'Gain re. psv (dB)')

    n = 5;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.3910,0.5930,0.2210,0.1600]);
    xlim(h(n),[-2.25 1])
    xticks(h(n),-2:0.5:1)
    xticklabels(h(n),{'','','','','','',''})
    ylim(h(n),[-20 40])
    yticks(h(n),-20:20:40)
    yticklabels(h(n),{'','','',''})

    n = 6;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.6660,0.5930,0.2210,0.1600]);
    xlim(h(n),[-2.25 1])
    xticks(h(n),-2:0.5:1)
    xticklabels(h(n),{'','','','','','',''})
    ylim(h(n),[-20 40])
    yticks(h(n),-20:20:40)
    yticklabels(h(n),{'','','',''})

    n = 7;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.1125,0.3960,0.2210,0.1600]);
    xlim(h(n),[-2.25 1])
    xticks(h(n),-2:0.5:1)
    xticklabels(h(n),{'-2','','-1','','0','','1'})
    ylim(h(n),[-1 1])
    yticks(h(n),-1:0.5:1)
    xlabel(h(n),'Frequency re. CF (oct.)') 
    ylabel(h(n),'Phase re. psv (deg.)')

    n = 8;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.3910,0.3960,0.2210,0.1600]);
    xlim(h(n),[-2.25 1])
    xticks(h(n),-2:0.5:1)
    xticklabels(h(n),{'-2','','-1','','0','','1'})
    ylim(h(n),[-1 1])
    yticks(h(n),-1:0.5:1)
    yticklabels(h(n),{'','','','',''})
    xlabel(h(n),'Frequency re. CF (oct.)')

    n = 9;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.6660,0.3960,0.2210,0.1600]);
    xlim(h(n),[-2.25 1])
    xticks(h(n),-2:0.5:1)
    xticklabels(h(n),{'-2','','-1','','0','','1'})
    ylim(h(n),[-1 1])
    yticks(h(n),-1:0.5:1)
    yticklabels(h(n),{'','','','',''})
    xlabel(h(n),'Frequency re. CF (oct.)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yy = yBM(loc,a:b)./d_ref(a:b); 
    plot(h(1),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color','k','LineWidth',1.5);

    yy = yDC(loc,a:b)./d_ref(a:b); 
    plot(h(1),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[1,0.2,0.2],'LineWidth',1.5);

    yy = yRL(loc,a:b)./d_ref(a:b); 
    plot(h(1),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.8,0.8],'LineWidth',1.5);

    yy = yBMp(loc,a:b)./d_ref(a:b); 
    plot(h(1),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color','k','LineWidth',1.5,'LineStyle','--');

    yy = yDCp(loc,a:b)./d_ref(a:b); 
    plot(h(1),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[1,0.2,0.2],'LineWidth',1.5,'LineStyle','--');

    yy = yRLp(loc,a:b)./d_ref(a:b);
    plot(h(1),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.8,0.8],'LineWidth',1.5,'LineStyle','--');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yy = xDC(loc,a:b)./d_ref(a:b); 
    plot(h(2),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[1,0.2,0.2],'LineWidth',1.5);

    yy = xRL(loc,a:b)./d_ref(a:b); 
    plot(h(2),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.8,0.8],'LineWidth',1.5);

    yy = xDCp(loc,a:b)./d_ref(a:b); 
    plot(h(2),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[1,0.2,0.2],'LineWidth',1.5,'LineStyle','--');

    yy = xRLp(loc,a:b)./d_ref(a:b);
    plot(h(2),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.8,0.8],'LineWidth',1.5,'LineStyle','--');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yy = zDC(loc,a:b)./d_ref(a:b); 
    plot(h(3),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[1,0.2,0.2],'LineWidth',1.5);

    yy = zRL(loc,a:b)./d_ref(a:b); 
    plot(h(3),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.8,0.8],'LineWidth',1.5);

    yy = zDCp(loc,a:b)./d_ref(a:b); 
    plot(h(3),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[1,0.2,0.2],'LineWidth',1.5,'LineStyle','--');

    yy = zRLp(loc,a:b)./d_ref(a:b);
    plot(h(3),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.8,0.8],'LineWidth',1.5,'LineStyle','--');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ref = yBMp(loc,a:b);
    yy = yBM(loc,a:b)./ref; 
    zero_ = round(mean(unwrap(angle(yy))/(2*pi)));
    plot(h(4),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.0,0.0],'LineWidth',1.5);
    plot(h(7),log2(freq(a:b)/freq(peak_idx)),unwrap(angle(yy))/(2*pi)-zero_,'Color',[0,0.0,0.0],'LineWidth',1.5);

    ref = yDCp(loc,a:b);          
    yy = yDC(loc,a:b)./ref; 
    zero_ = round(mean(unwrap(angle(yy))/(2*pi)));
    plot(h(4),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[1,0.2,0.2],'LineWidth',1.5);
    plot(h(7),log2(freq(a:b)/freq(peak_idx)),unwrap(angle(yy))/(2*pi)-zero_,'Color',[1,0.2,0.2],'LineWidth',1.5);

    ref = yRLp(loc,a:b);
    yy = yRL(loc,a:b)./ref; 
    zero_ = round(mean(unwrap(angle(yy))/(2*pi)));
    plot(h(4),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.8,0.8],'LineWidth',1.5);
    plot(h(7),log2(freq(a:b)/freq(peak_idx)),unwrap(angle(yy))/(2*pi)-zero_,'Color',[0,0.8,0.8],'LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ref = yDCp(loc,a:b);    
    yy = xDC(loc,a:b)./ref; 
    zero_ = round(mean(unwrap(angle(yy))/(2*pi)));
    pp = unwrap(angle(yy))/(2*pi)-zero_;
    if mean(pp) < -0.25
        pp = pp + 1;
    end
    plot(h(5),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[1,0.2,0.2],'LineWidth',1.5);
    plot(h(8),log2(freq(a:b)/freq(peak_idx)),pp,'Color',[1,0.2,0.2],'LineWidth',1.5);

    ref = xRLp(loc,a:b);    
    yy = xRL(loc,a:b)./ref;    
    zero_ = round(mean(unwrap(angle(yy))/(2*pi))); 
    pp = unwrap(angle(yy))/(2*pi)-zero_;
    if mean(pp) < -0.25
        pp = pp + 1;
    end
    plot(h(5),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.8,0.8],'LineWidth',1.5);
    plot(h(8),log2(freq(a:b)/freq(peak_idx)),pp,'Color',[0,0.8,0.8],'LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ref = zDCp(loc,a:b);
    yy = zDC(loc,a:b)./ref;
    zero_ = round(mean(unwrap(angle(yy))/(2*pi)));
    pp = unwrap(angle(yy))/(2*pi)-zero_;
    if mean(pp) < -0.25
        pp = pp + 1;
    end
    plot(h(6),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[1,0.2,0.2],'LineWidth',1.5);
    plot(h(9),log2(freq(a:b)/freq(peak_idx)),pp,'Color',[1,0.2,0.2],'LineWidth',1.5);

    ref = zRLp(loc,a:b);
    yy = zRL(loc,a:b)./ref;
    zero_ = round(mean(unwrap(angle(yy))/(2*pi)));
    pp = unwrap(angle(yy))/(2*pi)-zero_;
    if mean(pp) < -0.25
        pp = pp + 1;
    end
    plot(h(6),log2(freq(a:b)/freq(peak_idx)),20*log10(abs(yy)),'Color',[0,0.8,0.8],'LineWidth',1.5);
    plot(h(9),log2(freq(a:b)/freq(peak_idx)),pp,'Color',[0,0.8,0.8],'LineWidth',1.5);

    for i = 1:9%numel(h)
        naxis(h(i));
        xline(h(i),log2(freq(f_tail)/freq(peak_idx)),'Color','k','LineWidth',1.5);        
        xline(h(i),log2(freq(f_peak)/freq(peak_idx)),'Color','k','LineWidth',1.5);
    end

    text(h(1),log2(freq(f_tail)/freq(peak_idx))-0.3,55,'f1')
    text(h(1),log2(freq(f_peak)/freq(peak_idx))-0.3,55,'f2')
    text(h(2),log2(freq(f_tail)/freq(peak_idx))-0.3,55,'f1')
    text(h(2),log2(freq(f_peak)/freq(peak_idx))-0.3,55,'f2')
    text(h(3),log2(freq(f_tail)/freq(peak_idx))-0.3,55,'f1')
    text(h(3),log2(freq(f_peak)/freq(peak_idx))-0.3,55,'f2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pfact = 0.175;
    n = 10;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.030,0.140,0.2210,0.1600]);   

    l = plot_VC_section(h(n),MP,Nd,El,R(2),MP.xx(loc),f_tail,pfact,'-');
    text(h(n),0.4*l,0.2*l,'f1')
    h(n).Visible = 'off';

    n = 11;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.260,0.140,0.2210,0.1600]);   

    l = plot_VC_section(h(n),MP,Nd,El,R(2),MP.xx(loc),f_peak,pfact,'-');
    text(h(n),0.4*l,0.2*l,'f2')
    h(n).Visible = 'off';

    pfact = 0.16;
    n = 12;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.490,0.140,0.2210,0.1600]);   

    l = plot_VC_section(h(n),MP,Nd,El,R(1),MP.xx(loc),f_tail,pfact,'--');
    text(h(n),0.4*l,0.2*l,'f1')
    h(n).Visible = 'off';

    n = 13;
    h(n) = axes;
    hold(h(n),'on')
    set(h(n),'Position',[0.720,0.140,0.2210,0.1600]);   

    l = plot_VC_section(h(n),MP,Nd,El,R(1),MP.xx(loc),f_peak,pfact,'--');
    text(h(n),0.4*l,0.2*l,'f2')
    h(n).Visible = 'off';

end

function l = plot_VC_section(ax,MP,Nd,El,R,loc,nf,pfact,l_style)
% loc: location
% nf: # reference of the frequency

    NDOF = 6;

    nz = knnsearch(MP.xx(:),loc);
    nn = ((1:Nd.Nc)-1)*Nd.Nr + nz;

    xx = Nd.X(nz,:);
    yy = Nd.Y(nz,:);

    udof = MP.dof(2).udof;
    uu = reshape(R.Uf(udof,nf),NDOF,Nd.N);

    dx = uu(1,nn);
    dy = uu(2,nn);
    dz = uu(3,nn);
    
    ps = pfact*max(yy)/max(sqrt(abs(dx).^2 + abs(dy).^2));

    hold(ax,'on');
    n_sec = find(abs(Nd.Z + MP.loc*1e3 - loc*1e3) < 15);
    e_idx = find(ismember(El.Nd1,nn));

    nAP = find(strcmp(Nd.name,'AP'));
    nAM = find(strcmp(Nd.name,'AM'));
    nBB = find(strcmp(Nd.name,'BB'));
    nDD = find(strcmp(Nd.name,'DD'));
    nCC = find(strcmp(Nd.name,'CC'));
    nodes = [nAM(nz),nBB(nz),nDD(nz)];

    np = 21;
    phi = linspace(-pi,pi,np) + pi/np*rand(1,np);  
    icount = 0;
    for si = phi
        icount = icount + 1;
        for ie=1:El.N   
            if ismember(ie,e_idx)
            nds = [El.Nd1(ie), El.Nd2(ie)];               
                %undeformed shadow
                [R,T] = rotation_matrix(El,ie);
                d = uu(:,nds)*exp(1i*si);
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
                plot(ax,xx_l(1,:) + dd(1,:), xx_l(2,:) + dd(2,:),'Color',[0.6,0.6,0.6,0.2],'LineWidth', 1.5);                
            end
        end
    end
    
    np = 61;
    phi = linspace(0,2*pi,np) - angle(uu(2,nAM(nz)));  
    % cmap = hsv(length(phi));
    cmap = [0,  0,  0;
            0,0.8,0.8;
            1,0.2,0.2];
    icount = 0;
    for node = nodes
        % nodes
        icount = icount + 1;
        x = [Nd.X(node) + ps*real(uu(1,node)*exp(1i*phi)),NaN];
        y = [Nd.Y(node) + ps*real(uu(2,node)*exp(1i*phi)),NaN];

        plot(ax,x,y,'Color',cmap(icount,:),'LineWidth',2.5,'LineStyle',l_style);
        % patch(x,y,[phi,NaN],'EdgeColor','interp','LineWidth',3);
        % colormap(ax,cmap);
        % scatter(x,y,12,[cmap(icount,:)],'filled','o','Parent',ax);
    end
    xpadding = 0.2*(max(xx) - min(xx));
    ypadding = 0.2*(max(yy) - min(yy));
    axis(ax,'equal')
    xlim([min(xx) + xpadding, max(xx) - xpadding]); % should be adjusted manually for now
    ylim([min(yy) - ypadding, max(yy) + ypadding]); % should be adjusted manually for now

    l = max(xx) - min(xx);
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

    O33 = zeros(3);    
    R = [T, O33, O33, O33;  O33, T, O33, O33;  O33, O33, T, O33;   O33, O33, O33, T;];
end