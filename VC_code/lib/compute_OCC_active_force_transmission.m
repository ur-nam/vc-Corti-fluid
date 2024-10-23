function k = compute_OCC_active_force_transmission(loc,opt_plot)

% % xx = location where stiffnesses are computed [mm]
% % kOHC = axial stiffness of the OHC  [mN/m]
% % kOHB = stiffness of shear of the OHB  [mN/m]

    if ~exist('opt_plot','var')
        opt_plot = 0;
    end
%     set_Global_coefficients; % this is setting coefficients for parameter study
    lambdaF = -1;    % if it is zero, a point force. if negative uniformly distributed line force     
    MP.loc = loc;             % [mm]
    MP.length_BM = 600;   % [um]

    dim = 3;                    % dimension 2 or 3-D
    opt_lgradient = 0;          % gradient along the longitudinal direction
    opt_MP_set = 1;         % 0: the MP set before 2021.04

    MP.dZ = 10;
    MP.xx = MP.loc-0.5*(1e-3*MP.length_BM):(1e-3*MP.dZ):MP.loc+0.5*(1e-3*MP.length_BM);
    [Nd, El, MP] = model_3D(MP,dim,opt_lgradient, opt_MP_set);  

    nAP = find(strcmp(Nd.name,'AP'));
    nAM = find(strcmp(Nd.name,'AM'));
    nDD = find(strcmp(Nd.name,'DD'));
    nBB = find(strcmp(Nd.name,'BB'));
    eOHC = strcmpi(El.name,'OHC');
    eDCb = strcmpi(El.name,'DCb');
    eANK = strcmpi(El.name,'ANK');

    kOHC = 1e-3*El.YM(eOHC).*El.A(eOHC)./El.oL(eOHC);
    kDCb = 1e-3*El.YM(eDCb).*El.A(eDCb)./El.oL(eDCb);
    kOHB = 1e-3*El.YM(eANK).*El.A(eANK)./El.oL(eANK);    
    kOHB = 10*kOHB; % For this factor, refer to model_Moduli5d.m   

    idx = (MP.xx == MP.loc);

    % OHC force w BM fixed     
    Nd = set_boundary_conditions(MP, Nd, opt_MP_set, 1);              
    [rK,El] = createK(Nd, El); % reduced K matrix
    LU = createL(Nd, El);
    rK = rK + 1e9*LU.'*LU;
    K = createK_full(Nd,El); 
           
    fOHC = 2; fMET = 0;
    F_OHC = apply_active_force(lambdaF, Nd, El, MP, fOHC, fMET);
    b = F_OHC;
    a = rK\b;
    U_OHC = a(1:Nd.rdof);
  
    U_OHC_full = zeros(Nd.tdof,1);
    U_OHC_full(Nd.mapR2F(1:Nd.rdof)) = U_OHC;    

    if opt_plot
        figure; set(gcf,'Position',[30, 64, 445, 626]); clf;
        ax = axes;
    end

    Nd = update3D(Nd, U_OHC);     El = get_vector(El, Nd, 1);
    F = transpose(reshape(K*U_OHC_full,6,Nd.N)); % [dx,dy,dz,rx,ry,rz]
    dir = El.dir(eOHC,:);

    u_bot = dir(idx,:)*[Nd.dx(nDD(idx));Nd.dy(nDD(idx));Nd.dz(nDD(idx))];    
    u_top = dir(idx,:)*[Nd.dx(nBB(idx));Nd.dy(nBB(idx));Nd.dz(nBB(idx))];
    f_top = dir(idx,:)*F(nBB(idx),1:3).';
    f_bot = dir(idx,:)*F(nDD(idx),1:3).';

    ee = find(eDCb); ee = ee(idx);
    [ua,ub,fa,fb] = eval_disp_force(Nd,El,U_OHC_full,ee); 
    kDCj_a = 1e-3*(fa(2))/(ua(2)-ua(1));
    kDCj_b = 1e-3*(fb(2)-fb(1))/(ub(2)-ub(1));

%     F_AP = F(nAP(idx),1:3);
%     F_AM = F(nAM(idx),1:3);
    if opt_plot
%     plot_cross_section(ax,MP,Nd,El,F_OHC)
    plot_cross_section(subplot(2,1,1),MP,Nd,El,F);
    end

    kOCC_top = 1e-3*(f_top + kOHC*(u_bot-u_top))./u_top; % [nN/pN]*[pN]/[um] = [mN/m]
    kOCC_bot = 1e-3*(f_bot + kOHC*(u_top-u_bot))./u_bot; % in [mN/m]

    % OHC force w BM free
    Nd = set_boundary_conditions(MP, Nd, opt_MP_set, 0);
    [rK,El] = createK(Nd, El); % reduced K matrix
    LU = createL(Nd, El);
    rK = rK + 1e9*LU.'*LU;
    fOHC = 2; fMET = 0;    
    F_OHC = apply_active_force(lambdaF, Nd, El, MP, fOHC, fMET);
    b = F_OHC;
    a = rK\b;
    U_OHC = a(1:Nd.rdof);

    Nd = update3D(Nd, U_OHC);     El = get_vector(El, Nd, 1);
    xOHC = abs(El.L(eOHC)-El.oL(eOHC));
    kOCC_OHC = 1e-3*(1000)./xOHC - kOHC; % in [mN/m] % force is hardcoded    

    if opt_plot
%     plot_cross_section(ax,MP,Nd,El,F_OHC)
    plot_cross_section(subplot(2,1,2),MP,Nd,El,F_OHC);
    end
    k.OCC_OHC = kOCC_OHC(idx);
    k.OCC_top = kOCC_top(idx);
    k.OCC_bot = kOCC_bot(idx);

    k.rOCC_top = kOCC_top(idx)./kOHC(idx);
    k.rOCC_bot = kOCC_bot(idx)./kOHC(idx);
    k.rOCC_OHC = kOCC_OHC(idx)./kOHC(idx);

    k.DCj_a = kDCj_a;
    k.DCj_b = kDCj_b;
    
    k.OHC = kOHC(idx);
    k.OHB = kOHB(idx);

end

function lu = createL(Nd, El)

    NDOF = 6;
    udof = Nd.tdof;
    ldof = 3*0.5*numel(Nd.Master_Slave);

    partsR = zeros(numel(Nd.Master_Slave)*3,4);
    cnt = 0;

    for ni = 1:length(Nd.Master_Slave)
        m = Nd.Master_Slave(ni,1);
        s = Nd.Master_Slave(ni,2);
        m_dof = (m-1)*NDOF + (1:3);
        s_dof = (s-1)*NDOF + (1:3);
        l_dof = (ni-1)*3 + (1:3);

        rdof = l_dof;
        cdof = [m_dof,s_dof];
        m = [1,0,0,-1,0,0;
             0,1,0,0,-1,0;
             0,0,1,0,0,-1;];
        [partsR,cnt] = fillR(rdof,cdof,m,partsR,cnt);
    end

    lu = sparse(partsR(:,1),partsR(:,2),partsR(:,3),ldof,udof);
    lu = lu(:,Nd.BC);

end

function K = createK_full(Nd,El)

    NDOF = 6;
    partsK = zeros(144*El.N,3);
    cnt = 0; % counter variable

    for i = 1:El.N
    
        nd1 = El.Nd1(i);
        nd2 = El.Nd2(i);
        k = El.mat(i).k;
        % for the beam elements
            if El.type(i) < 1

                rdof = [(nd1-1)*NDOF + (1:6),(nd2-1)*NDOF + (1:6)];
                
                % This function fills in the 3-column input for the sparse
                % function.
                [partsK,cnt] = fillR(rdof,rdof,k,partsK,cnt);
                
            else
                rdof = [(nd1-1)*NDOF + (1:3),(nd2-1)*NDOF + (1:3)];                
                [partsK,cnt] = fillR(rdof,rdof,k,partsK,cnt);
                
            end
    end

    % Remove excess zero elements resulting from 144 DoF assumption.
    zr = partsK(:,1)==0;
    partsK(zr,:) = [];
    % Use sparse function to create the matrix
    K = sparse(partsK(:,1),partsK(:,2),partsK(:,3),Nd.tdof,Nd.tdof);
end

function plot_cross_section(ax,MP,Nd,El,F)

    pfact = 10;

    nz = 0.5*MP.NZ + 1;
    xx = Nd.X(nz,:);
    yy = Nd.Y(nz,:);
    zz = Nd.Z(nz,:);

    dx = Nd.dx(nz,:);
    dy = Nd.dy(nz,:);
    dz = Nd.dz(nz,:);
    rx = Nd.rx(nz,:);
    ry = Nd.ry(nz,:);
    rz = Nd.rz(nz,:);

    uu = transpose([Nd.dx(:),Nd.dy(:),Nd.dz(:),Nd.rx(:),Nd.ry(:),Nd.rz(:)]);

    ps = pfact/max([abs(dx),abs(dy)]);

    nn = ((1:Nd.Nc)-1)*Nd.Nr + nz;
    hold(ax,'on');
    for ie=1:El.N        
        nds = [El.Nd1(ie), El.Nd2(ie)];       
        if ismember(nds,nn)
            plot(ax,Nd.X(nds), Nd.Y(nds), 'Color',[0.8 0.8 0.8],'LineWidth', 1.0);
            if El.type(ie) < 1
                [R,T] = rotation_matrix(El,ie);
                d = uu(:,nds);
                coord = [Nd.X(nds);Nd.Y(nds)];
                d = R*d(:); % move displacement values to local frame
                L = El.L(ie);
                nl = 21;
                x = linspace(0,El.L(ie),nl);
                xx_l = zeros(2,nl);
                dy_l = zeros(1,nl);
                dx_l = zeros(1,nl);
                O = zeros(1,nl);
                for ii = 1:numel(x)
                    % this is in x-y plane
                    N = [(L-x(ii))/L, x(ii)/L]; %[u1,u2]
                    xx_l(:,ii) = coord*N.';
                    dx_l(ii) = (N*d([1;7]));
                    N = [1-3*x(ii)^2/L^2+2*x(ii)^3/L^3, x(ii)-2*x(ii)^2/L+x(ii)^3/L^2,...
                         3*x(ii)^2/L^2-2*x(ii)^3/L^3,-x(ii)^2/L+x(ii)^3/L^2]; %[v1,t1,v2,t2] local
                    dy_l(ii) = (N*d([2;6;8;12]));
                end
                dd = T.'*ps*[dx_l;dy_l;O];
                plot(ax,xx_l(1,:) + dd(1,:), xx_l(2,:) + dd(2,:), 'Color',[0 0 0],'LineWidth', 2.5);
            else
                plot(ax,Nd.X(nds) + ps*Nd.dx(nds), Nd.Y(nds) + ps*Nd.dy(nds), 'Color',[0 0 0],'LineWidth', 2.5);
            end
        end
    end

    plot(ax,xx + ps*dx, yy + ps*dy, 'o','Color','none','MarkerFaceCol',[0, 0.4470, 0.7410],'MarkerSize',4);

    tmp = zeros(Nd.tdof,1);
    if numel(F) == Nd.rdof
        tmp(Nd.mapR2F(1:Nd.rdof)) = F;
        F = 3e-2*reshape(tmp,6,Nd.N);
    else
        F = 3e-2*transpose(F);
    end
    fx = F(1,nn);
    fy = F(2,nn);

    quiver(ax,xx + ps*dx, yy + ps*dy, fx, fy,'LineWidth',2,...,
        'Autoscale','off','Color','r');
    daspect(ax,[1 1 1]);
    maxF = max(abs(F));

    padding = 0.1*(max(xx) - min(xx));
    xlim([min(xx) - padding, max(xx) + padding]); % should be adjusted manually for now
    ylim([min(yy) - padding, max(yy) + padding]); % should be adjusted manually for now

end

function [xx, kTM] = compute_TM_axial_stiffness(Nd, El)

        xx = (min(Nd.Z(:,1)):10:max(Nd.Z(:,1)))';
        n = length(xx);
        kTM = zeros(size(xx));

        for ii = 1:n

            xi = xx(ii);

            eTMx = find(strcmpi(El.name, 'TMx1') | strcmpi(El.name, 'TMx2'));
            eTMx = eTMx(Nd.Z(El.Nd1(eTMx))< xi + 0 & Nd.Z(El.Nd1(eTMx))>= xi-5);    

            ce = El.oL(eTMx)./(El.YM(eTMx).*El.A(eTMx));
            kTM(ii) = 1/sum(ce);
        end

end

function [xx, kTM] = compute_TM_bending_stiffness(Nd, El)

    xx = (min(Nd.Z(:,1)):10:max(Nd.Z(:,1)))';
    n = length(xx);
    kTM = zeros(size(xx));

    for ii = 1:n
        
        xi = xx(ii);

        eTMx1 = find(strcmpi(El.name, 'TMx1'));   % attachment section
        eTMx2 = find(strcmpi(El.name, 'TMx2'));   % body section

        eTMx1 = eTMx1(Nd.Z(El.Nd1(eTMx1))< xi+0 & Nd.Z(El.Nd1(eTMx1))>=xi-5);
        eTMx2 = eTMx2(Nd.Z(El.Nd1(eTMx2))< xi+0 & Nd.Z(El.Nd1(eTMx2))>=xi-5);

        L1 = sum(El.oL(eTMx1));
        L2 = sum(El.oL(eTMx2));

        Iz1 = mean(El.Iz(eTMx1));
        Iz2 = mean(El.Iz(eTMx2));

        E1 = mean(El.YM(eTMx1));
        E2 = mean(El.YM(eTMx2));


        cTM = (L1^2)/(E1*Iz1)*(5/6*L1 + 3/2*L2) + (L2^3)/(3*E2*Iz2);

        kTM(ii) = 1/cTM;
    end

end
   
    
function [F, MP] = apply_BM_force(lambdaF, Nd, MP, maxF, opt)
    
    fdir = 2;
    
    fnode = find(strcmpi(Nd.name,'AM'));
    NE = length(fnode);
    ZZ = (1:NE)-(NE+1)/2;    

    if lambdaF >0
        sz = lambdaF/MP.dZ;
        ff = 1/sqrt(2*pi)/sz*exp(-(ZZ).^2/2/sz^2);
        ff = ff/sum(ff);
    elseif lambdaF == 0 % point force              
        ff = zeros(size(ZZ));
        ff(ZZ==0) = 0.5;
        ff(abs(ZZ)==1) = 0.25;
    elseif lambdaF < 0 % equally distributed force along the mid-BM,
        ff = ones(NE,1);
        
        if lambdaF == -999  % uniformly distributed pressure
            fnode = find(Nd.Y < 1e-3*MP.dZ & ~strcmpi(Nd.name,'AX') & ~strcmpi(Nd.name,'A0')); % entire BM
            NE = length(fnode);
            ff = 1e-3*ones(NE,1).*Nd.pArea(fnode);   % 1e-3 [pN/um^2] = 1e-3 [Pa] = 1 [mPa]
        end
        
    else  
        msg = 'Error in ''apply_BM_force''';
        fprintf(1,['\n' msg '\n']);        
        ff = zeros(size(ZZ));
    end
    
        
    MP.fBM = maxF*ff;
    
    if exist('opt','var')
        if opt == 0 % pressure distributed evenly along radial section
            fnode = find(Nd.Y==0 & Nd.X~=0 & ~strcmp(Nd.name,'AX'));  % entire BM            
            NE = length(fnode);
            ff = zeros(1,NE);
            sz = lambdaF;
            for ii=1:NE
                zi = Nd.Z(fnode(ii));
                ff(ii) = 1/sqrt(2*pi)/sz*exp(-(zi).^2/2/sz^2);                
            end
%             ff = ff/sum(ff);            
        end
    end

    F = 0;
    
    if maxF ~= 0
        maxff=max(ff);
        for ke = 1:length(fnode)
            if abs(ff(ke)/maxff) > 0.001
                F = F + maxF*ff(ke)*make_unit_force2(Nd, fnode(ke),fdir);
            end
        end
    end
    

end

    
function [F, MP] = apply_active_force(lambdaF, Nd, El, MP, fOHC, fMET)

       
        sz = lambdaF/MP.dZ;
                        
        ie_ohc = find( strcmp('OHC',El.name) ==1);
        NE = length(ie_ohc);
        ZZ = (1:NE)-(NE+1)/2;
        
        if lambdaF >0            
            ff = 1/sqrt(2*pi)/sz*exp(-(ZZ).^2/2/sz^2);
            ff = ff/max(ff);
        elseif lambdaF == 0 % point force              
            ff = zeros(size(ZZ));
            ff(ZZ==0) = 1.0;
        elseif lambdaF < 0 % equally distributed force along the length,
            ff = ones(NE,1);
        else  
            msg = 'Error in ''apply_active_force''';
            fprintf(1,['\n' msg '\n']);        
            ff = zeros(size(ZZ));
        end        

        if fOHC ~= 0            % Force by OHC prestine
            
            nd1 = El.Nd1(ie_ohc); %dd
            nd2 = El.Nd2(ie_ohc); %bb
                                      
            F = 0;            
            %qOHC = fOHC*1e3;   % 1000 pN per 10 um section
            qOHC = 1*1e3;   % 1000 pN per 10 um section
            for ke = 1:length(nd1)  
                ie = ie_ohc(ke);
                fdir = El.dir(ie,:);
                if (fOHC == 2) || (fOHC == 1)
                    F = F - qOHC*ff(ke)*make_unit_force2(Nd,nd2(ke),fdir);
                end
                if (fOHC == 2) || (fOHC == -1)
                    F = F + qOHC*ff(ke)*make_unit_force2(Nd,nd1(ke),fdir);
                end                
            end
            
        else
            qOHC = 0;            
            F = 0;
        end
        
        MP.fOHC = qOHC*ff;     % applied somatic force per HC  
        
        if fMET ~= 0            % Force by MET channel
            
            ie_ohc = find(strcmp('OHB',El.name)==1);
            nd1 = El.Nd1(ie_ohc);
            nd2 = El.Nd2(ie_ohc);
            je_ohc = find(strcmp('ANK',El.name)==1);

            
            qMET = fMET*1e3;
            
            
            for ke=1:NE
                fdir = El.dir(je_ohc(ke),:);                
                F = F - qMET*ff(ke)*make_unit_force2(Nd,nd1(ke),fdir);
                F = F + qMET*ff(ke)*make_unit_force2(Nd,nd2(ke),fdir);                
            end
        else
            qMET = 0;
        end 
        
        MP.fMET = -qMET*ff/3;     % applied hair bundle shear force per HC
end

function Nd = set_boundary_conditions(MP, Nd, opt_MP_set, opt_BM_fix)

    NDOF = 6;
    if ~exist('opt_MP_set','var') %  a new set of geometry and material props updated on Apr 2021
        opt_MP_set = 1;         % default is using the new set. 0 = the set before 2021
    end

    Nd.BC = ones(Nd.N,NDOF);

    nd_SLam = find(strcmp(Nd.name,'A0'));       % Spiral laminar at IPC foot (medial edge)
    nd_SLig = find(strcmp(Nd.name,'AX'));       % Spiral ligament (lateral edge)
    nd_TLam = find(strcmp(Nd.name,'E0'));       % Spiral lamina at TM attach't    
    nd_BM = find(strncmp(Nd.name,'A',1));       % BM nodes

    nd_IPCToe = find(strcmp(Nd.name,'A1'));     % IPC toe
%     nd_IPCToe = find(strcmp(Nd.name,'A2'));     % IPC toe

    [nr,nc] = size(Nd.X);    
%     nd_bend = (1:nr:nr*nc);
    locz = -0.5*(MP.length_BM); % [um]
    nd_bend = find(Nd.Z<0.4*MP.dZ + locz & Nd.Z>=-0.5*MP.dZ + locz);
    locz = 0.5*(MP.length_BM); % [um]
    nd_aend = find(Nd.Z<0.4*MP.dZ + locz & Nd.Z>=-0.5*MP.dZ + locz);
%     nd_aend = (nr:nr:nr*nc);
    %nd_aend = [(nr-1:nr:(nr-1)*nc), (nr:nr:nr*nc)];
    nd_ends = [nd_bend; nd_aend];
%     nd_ends = [nd_aend];
    if opt_MP_set == 0    
        Nd = set_BC3(Nd, nd_SLam, 'dAll');
    else
        Nd = set_BC3(Nd, nd_SLam, 'dAll');        
        Nd = set_BC3(Nd, nd_IPCToe, 'dAll');
    end
    Nd = set_BC3(Nd, nd_SLig, 'All');
    
    Nd = set_BC3(Nd, nd_TLam, 'All');

    Nd = set_BC3(Nd, nd_ends, 'All'); 
%     Nd = set_BC3(Nd, nd_ends, 'dz');    
%     Nd = set_BC3(Nd, nd_ends, 'dx');
%     Nd = set_BC3(Nd, nd_ends, 'ry');
%     Nd = set_BC3(Nd, nd_ends, 'rz');  
    if opt_BM_fix, Nd = set_BC3(Nd, nd_BM, 'All'); end

    Nd = find_mapF2R2(Nd);
    
    % MSA: this was adjusted in the end to make using Nd.BC easier 
    Nd.BC = reshape(transpose(logical(Nd.BC)),Nd.tdof,1);
    
end



function Nd = set_BC3(Nd, cnode, cdof)

% %      Nd = set_boundary_condition(Nd, cnode, cdof)
% % 
% %      Description: assign boundary conditions
% % 
% %      Nd: node structure, cnode: node id to constrain, 
% %      cdof: DOF to confine (1=dx, 4=rx, ...)
% %     
% %      First written: 2009.03.27
% %      Last modified: 
% %      by Jong-Hoon Nam
% % 

    if ischar(cnode)
        if strcmpi(cnode(1),'A')
            cnode = 1:Nd.N;
        end
    end
    
    if ~isempty(cnode)
    
        switch (upper(cdof(1)))
            case 'A'
                Nd.BC(cnode,:) = 0;
            case 'D'
                if strcmpi(cdof(2),'X')
                    Nd.BC(cnode,1) = 0;
                elseif strcmpi(cdof(2),'Y')
                    Nd.BC(cnode,2) = 0;
                elseif strcmpi(cdof(2),'Z')
                    Nd.BC(cnode,3) = 0;
                elseif strcmpi(cdof(2),'A')
                    Nd.BC(cnode,1:3) = 0;
                else
                    printf(1,'\n Warning: Unknown option for displ constraint \n');
                end
            case 'R'
                if strcmpi(cdof(2),'X')
                    Nd.BC(cnode,4) = 0;
                elseif strcmpi(cdof(2),'Y')
                    Nd.BC(cnode,5) = 0;
                elseif strcmpi(cdof(2),'Z')
                    Nd.BC(cnode,6) = 0;
                elseif strcmpi(cdof(2),'A')
                    Nd.BC(cnode,4:6) = 0;
                else
                    printf(1,'\n Warning: Unknown option for rotational constraint \n');
                end            
            otherwise
                fprintf(1,'\n Warning: check if the boundary condition was properly assigned. \n');
        end        
    end
end

function [ua,ub,fa,fb] = eval_disp_force(Nd,El,Uf,ee)

    NDOF = 6;
    nd1 = El.Nd1(ee); ndof1 = (nd1-1)*NDOF + (1:6);
    nd2 = El.Nd2(ee); ndof2 = (nd2-1)*NDOF + (1:6);
    
    dir = El.dir(ee,:);
    nd = [nd1,nd2];
    ua = dir*[Nd.dx(nd);Nd.dy(nd);Nd.dz(nd)];
    ub = [-dir(2),dir(1),dir(3)]*[Nd.dx(nd);Nd.dy(nd);Nd.dz(nd)];

    f_ee = reshape(El.mat(ee).k*Uf([ndof1,ndof2]),NDOF,2);
    fa = dir*f_ee(1:3,:);
    fb = [-dir(2),dir(1),dir(3)]*f_ee(1:3,:);

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