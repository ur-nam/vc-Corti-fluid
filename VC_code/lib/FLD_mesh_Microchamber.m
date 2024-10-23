function FLD = FLD_mesh_Microchamber(MP,H,FLD)
    H_05 = round(H/2);
    H_OCC = 25;
    L = MP.length_BM;
    Lo = 1000;
    Ho = 80;
    Xo = 2500;
    Ls = 1000;
    Hs = 200;
    Xs = 2500;
    Heli = 100;
    Lt = L+Heli;
    XC1 = L+Heli-H_05;
    XC2 = L;
    botH = 1000;
    topH = 500;
%     Slit = ['drectangle(p,0,1000,-200,0)'];
%     Topc = ['drectangle(p,-1000,2000,0,1000)'];
%     CP = ['drectangle(p,0,1000,-10,10)'];
%     Botc = ['drectangle(p,-2000,3000,-1200,-200)'];
    Gap = ['drectangle(p,',num2str(Xo),',',num2str(Xo+Lo),',',num2str(H/2),',',num2str(H/2+Ho),')'];
    Topc = ['drectangle(p,',num2str(Xo-1000),',',num2str(Xo+Lo+1000),',',num2str(H/2+Ho),',',num2str(H/2+Ho+topH),')'];
    CP = ['drectangle(p,0,',num2str(L),',',num2str(-H_OCC),',',num2str(H_OCC),')'];
    CD = ['drectangle(p,0,',num2str(L),',',num2str(-H/2),',',num2str(H/2),')'];
    Botc = ['drectangle(p,',num2str(Xs-1000),',',num2str(Xs+Ls+1000),',',num2str(-H/2-Hs-botH),',',num2str(-H/2-Hs),')'];
    Slit = ['drectangle(p,',num2str(Xs),',',num2str(Xs+Ls),',',num2str(-H/2-Hs),',',num2str(-H/2),')'];
    
%     a = ['(drectangle(p,0,',num2str(L),',',num2str(-H_05),',',num2str(H_05),'))'];
%     funinline = ['(drectangle(p,0,',num2str(L),',',num2str(-H_05),',',...
%         num2str(H_05),'))'];
%     funinline = ['ddiff(drectangle(p,0,',num2str(Lt),',',num2str(-H_05),',',...
%         num2str(H_05),'),','drectangle(p,',num2str(Heli),',',num2str(L),','...
%         ,num2str(-H_OCC),',',num2str(H_OCC),'))'];
% funinline = ['ddiff(dunion(drectangle(p,0,',num2str(XC1),',',num2str(-H_05),',',...
%         num2str(H_05),'),','dcircle(p,',num2str(XC1),',0,',num2str(H_05),')),','dunion(drectangle(p,0,',num2str(L),','...
%         ,num2str(-H_OCC),',',num2str(H_OCC),'),','dcircle(p,',num2str(XC2),',0,',num2str(H_OCC),')))'];
% funinline = ['ddiff(dunion(dunion(',Slit,',',Topc,'),',Botc,'),',CP,')'];

funinline = ['ddiff(dunion(dunion(dunion(dunion(',Topc,',',Gap,'),',CD,'),',Botc,'),',Slit,'),',CP,')']; 
% funinline = ['dunion(',Topc,',',Slit,')'];
    fd=inline(funinline,'p');
%     nz = round(L/10)+1;
nz = round(L/10)+1;
% fix nodes at BM and TM    
    for ii = 1:nz
        pfix(ii,1)= (ii-1)*10;
        pfix(ii,2)= H_OCC;
        pfix(ii+nz,1) = (ii-1)*10;
        pfix(ii+nz,2) = -H_OCC;
    end;
% fix nodes below BM and above TM
    for ii = 1:nz
        pfix(ii+2*nz,1)= (ii-1)*10;
        pfix(ii+2*nz,2)= H_OCC+10;
        pfix(ii+3*nz,1) = (ii-1)*10;
        pfix(ii+3*nz,2) = -H_OCC-10;
    end;
%     pfix(4*nz+1,:) = [-2500,-500];
%     pfix(4*nz+2,:) = [-2500,-1500];
%     pfix(4*nz+3,:) = [3000,-500];
%     pfix(4*nz+4,:) = [3000,-1500];


%     pfix(4*nz+1,:) = [-1000,1000];
%     pfix(4*nz+2,:) = [2000,1000];
%     pfix(4*nz+3,:) = [-2000,-1200];
%     pfix(4*nz+4,:) = [3000,-1200];
%     pfix(4*nz+5,:) = [-1000,0];
%     pfix(4*nz+6,:) = [2000,0];
%     pfix(4*nz+7,:) = [-2000,-200];
%     pfix(4*nz+8,:) = [3000,-200];
%     pfix(4*nz+9,:) = [0,-200];
%     pfix(4*nz+10,:) = [1000,-200];
    
    pfix(4*nz+1,:) = [Xo-1000,H/2+Ho+topH];
    pfix(4*nz+2,:) = [Xo+Lo+1000,H/2+Ho+topH];
    pfix(4*nz+3,:) = [0,-H/2];
    pfix(4*nz+4,:) = [L,-H/2];
    pfix(4*nz+5,:) = [Xo-1000,H/2+Ho];
    pfix(4*nz+6,:) = [Xo+Lo+1000,H/2+Ho];
    pfix(4*nz+7,:) = [Xo,H/2];
    pfix(4*nz+8,:) = [Xo+Lo,H/2];
    pfix(4*nz+9,:) = [0,H/2];
    pfix(4*nz+10,:) = [L,H/2];
    pfix(4*nz+11,:) = [Xs-1000,-H/2-Hs-botH];
    pfix(4*nz+12,:) = [Xs+Ls+1000,-H/2-Hs-botH];
    pfix(4*nz+13,:) = [Xs-1000,-H/2-Hs];
    pfix(4*nz+14,:) = [Xs+Ls+1000,-H/2-Hs];
    pfix(4*nz+15,:) = [Xs,-H/2];
    pfix(4*nz+16,:) = [Xs+Ls,-H/2];
    pfix(4*nz+17,:) = [Xs,-H/2-Hs];
    pfix(4*nz+18,:) = [Xs+Ls,-H/2-Hs];
    pfix(4*nz+19,:) = [Xo,H/2+Ho];
    pfix(4*nz+20,:) = [Xo+Lo,H/2+Ho];
    
%     hinline = 'min(abs(p(:,2).^1.3)/100+(abs(p(:,1)-500).^2.2)/250000+2.5,50)';
% hinline = ['min(abs(p(:,2).^1.3)/100+(abs(p(:,1)-',num2str(round(L/2)),').^2.2)/250000+2.5,50)'];
hinline = ['min(abs(p(:,2).^1.4)/100+2.5,30)'];

% hinline = '10';
    fh = inline(hinline,'p');
    h0 = 8;
    [FLD.p,FLD.t] = distmesh2d(fd,fh,h0,[min(0,Xs-1000),-H/2-Hs-botH;max(L,Xo+Lo+1000),H/2+Ho+topH],pfix);
    FLD.h0 = h0;
    FLD.HOC = H_OCC;
    FLD.Lo = Lo;
    FLD.Xo = Xo;
    FLD.Ho = Ho;
    FLD.Ls = Ls;
    FLD.Xs = Xs;
    FLD.Hs = Hs;
    FLD.L = L;
    FLD.H = H;
    FLD.Lt = Lt;
    FLD.Heli = Heli;


end