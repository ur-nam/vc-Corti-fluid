function plot_output(h,MP,fMP,Nd,El,AMat,Uf,po0,V0)

    fig = h(6:7);

    section = [201,601,1001];

    odof = MP.dof(5).odof;
    po_f = Uf(odof,:);
    po = inverse_fourier_transform(po_f);
    udof = MP.dof(2).udof;
    uu = inverse_fourier_transform(Uf(udof,:));
    edof = MP.dof(3).edof;
    ee = inverse_fourier_transform(Uf(edof,:)); 
   
    plot_HB_disp(fig(1),MP,fMP,Nd,El,uu,section)

    plot_po(fig(2),MP,fMP,po,section)

%     plot_Vm(fig(2),MP,fMP,ee,section)

%     plot_fOHC(MP,fMP,Nd,El,-AMat.Aue,ee,V0,section)
% 
%     plot_fMET(MP,fMP,Nd,El,-AMat.Auo,po,po0,section)

end

function plot_fMET(MP,fMP,Nd,El,Auo,yy,yy0,nn)

    ndof = 6;
    idx = reshape(transpose((1:Nd.tdof)),ndof,Nd.N);
    nOHB = find(strcmpi('OHB',El.name));
    nANK = find(strcmpi('ANK',El.name));
    fdir = El.dir(nANK,:);
%         nd1OHB = El.Nd1(nOHB);
    nd2OHB = El.Nd2(nOHB);
    f = zeros(MP.dof(2).n,fMP.L); % full vector
    f(Nd.BC,:) = Auo*(yy-yy0); % update using reduced matrix
    dir = 2;
    fMET = f(idx(dir,nd2OHB),:)./fdir(:,dir); % you can project back using any of the three components
    loc = MP.xx(nn);
    figure('Name','MET force');
    set(gcf,'Position',[1400,50,400,400])
    for i = 1:length(nn)
        subplot(length(nn),1,i)
        plot(fMP.tt,fMET(nn(i),:));
        title([num2str(loc(i)),'mm'])
    end   
end

function plot_fOHC(MP,fMP,Nd,El,Aue,yy,yy0,nn)

    ndof = 6;
    idx = reshape(transpose((1:Nd.tdof)),ndof,Nd.N);
    nOHC = find(strcmpi('OHC',El.name));
    fdir = El.dir(nOHC,:);
%         nd1OHC = El.Nd1(nOHC);
    nd2OHC = El.Nd2(nOHC);
    f = zeros(MP.dof(2).n,fMP.L); % full vector
    f(Nd.BC,:) = Aue*(yy-yy0); % update using reduced matrix
    dir = 2;
    fOHC = f(idx(dir,nd2OHC),:)./fdir(:,dir); % you can project back using any of the three components

    loc = MP.xx(nn);
    figure('Name','OHC force');
    set(gcf,'Position',[1400,550,400,400])
    for i = 1:length(nn)
        subplot(length(nn),1,i)
        plot(fMP.tt,fOHC(nn(i),:));
        title([num2str(loc(i)),'mm'])
    end   
end

function plot_Vm(fig,MP,fMP,yy,nn)

    nv = 5;

    V2 = (nn-1)*nv + 2;
    V3 = (nn-1)*nv + 3;

    Vm = yy(V2,:) - yy(V3,:);

    loc = MP.xx(nn);
    figure(fig); clf;
    set(fig,'Position',[950,550,400,400])
    for i = 1:length(nn)
        subplot(length(nn),1,i)
        plot(fMP.tt,Vm(i,:));
        title([num2str(loc(i)),'mm'])
    end      

end

function plot_po(fig,MP,fMP,yy,nn)

    loc = MP.xx(nn);
    figure(fig); clf;
    set(fig,'Position',[500,550,400,400])
    for i = 1:length(nn)
        subplot(length(nn),1,i)
        plot(fMP.tt,yy(nn(i),:));
        title([num2str(loc(i)),'mm'])
    end    

end

function plot_HB_disp(fig,MP,fMP,Nd,El,yy,nn)

    nOHB = find(strcmpi('OHB',El.name));
    nANK = find(strcmpi('ANK',El.name));
    nd2OHB = El.Nd2(nOHB);
    nd1OHB = El.Nd1(nOHB);
    nd1ANK = El.Nd1(nANK); 
    nd2ANK = El.Nd2(nANK);

    HBdir = 1;
    ndof = 6;

    dx2 = yy((nd2OHB - 1)*ndof + HBdir,:);
    dx1 = yy((nd1OHB - 1)*ndof + HBdir,:);
    dxOHB = 1e3*(dx2 - dx1);

    dx2 = yy((nd2ANK - 1)*ndof + HBdir,:);
    dx1 = yy((nd1ANK - 1)*ndof + HBdir,:);
    dxANK = 1e3*(dx2 - dx1);

    loc = MP.xx(nn);
    figure(fig); clf;
    set(fig,'Position',[50,550,400,400])
    for i = 1:length(nn)
        subplot(length(nn),1,i)
        plot(fMP.tt,dxOHB(nn(i),:));
        hold on
        plot(fMP.tt,dxANK(nn(i),:));
        legend({'OHB','ANK'});
        title([num2str(loc(i)),'mm'])
    end
        
end