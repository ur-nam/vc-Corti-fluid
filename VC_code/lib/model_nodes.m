function [Nd, MP] = model_nodes(MP)
% # function [Nd, MP] = model_nodes(MP)
% #
% # First written: 2009.06.xx
% # Last modified: 2013.03.29
% #
% # Written by: Jong-Hoon Nam at University of Rochester
% #
% # Assign nodal coordinates
% #   
        
    qMP = consider_geometrical_lgradient(MP);
    p = define_key_coordinates(MP, qMP);
    
    Nd = assign_nodal_coordinates(p, MP);
    
    Nd = initialize_nodes(Nd);
end

    
    
function p = define_key_coordinates(MP, qMP)

    if isfield(MP,'nz')
        nz = MP.nz;
    else
        nz = round(MP.length_BM/MP.dZ) + 1;
    end
          
    p.A0 = zeros(nz,2);                         % root of IPC
% %     p.A1 = [1/3*qMP.width_BM, zeros(nz,1)];     % root of OPC
% %     p.A2 = [1/2*qMP.width_BM, zeros(nz,1)];  % root of DC


    p.A1 = [qMP.OPC_root_location.*qMP.width_BM, zeros(nz,1)];     % root of OPC
    p.A2 = [qMP.DC_root_location.*qMP.width_BM, zeros(nz,1)];  % root of DC
    p.A3 = [qMP.width_BM,zeros(nz,1)];      % BM attachment spiral ligament    
    
    p.CC = [ qMP.width_TC, qMP.height_TC ];      % Tip of TC (tunnel of Corti)
    
    q0 = qMP.angle_RC*pi/180;
    p.BB = zeros(nz,2);
    for iz=1:nz
        q1 = q0(iz);
        R1 = [cos(q1), -sin(q1); sin(q1), cos(q1)];
        bb = R1*([qMP.width_BM(iz)/8;0]) + p.CC(iz,:)';   % root of OHC bundle or apex of OHC
               % Where it is assumed that the distance btw the tip of PC triangle and the 2nd row OHC bundle is 1/10 of the BM width
        pp = R1*([qMP.width_BM(iz)/8+qMP.diam_DCpr(iz);0]) + p.CC(iz,:)';   % apex of DC process
%         pp = R1*([qMP.width_BM(iz)/6+qMP.diam_DCpr(iz);0]) + p.CC(iz,:)';   % apex of DC process
        p.BB(iz,:) = bb';
        p.PP(iz,:) = pp';
    end    
 
    p.DD = zeros(nz,2);
    for iz=1:nz
        
        BB=p.BB(iz,:); 
        A2=p.A2(iz,:);        
        
        alpha = qMP.DD_alpha(iz);               
        theta = qMP.DD_theta(iz);
        
        DD = alpha*BB + (1-alpha)*A2;
        
        
        q1 = theta*pi/180;
        R1 = [cos(q1), -sin(q1); sin(q1), cos(q1)];
        

        p.DD(iz,:) = R1*(DD'-BB') + BB';   % joint btw OHC and DC
        
    end
    
    p.E1 = zeros(nz,2);
    p.FF = zeros(nz,2);
    p.E6 = zeros(nz,2);
    for iz=1:nz
        angle_RC = qMP.angle_RC(iz);
        q1 = (angle_RC-1)*pi/180;     R1 = [cos(q1), -sin(q1); sin(q1), cos(q1)];
        E1 = R1*[0; qMP.height_HB(iz)] + p.BB(iz,:)'; % TM - HB attachment point
        p.E1(iz,:) = E1';

        q1 = (angle_RC + 4)*pi/180;     R1 = [cos(q1), -sin(q1); sin(q1), cos(q1)];
        FF = R1*[0; qMP.height_HB(iz)] + p.BB(iz,:)'; % TM - HB attachment anchor point
        p.FF(iz,:) = FF';

        p.E6(iz,:) = [E1(1)-qMP.width_TM(iz), qMP.height_TMSLa(iz)];
    end
    
    p.E2 = 5/6*p.E1 + 1/6*p.E6;
    p.E3 = 4/6*p.E1 + 2/6*p.E6;
    p.E4 = 3/6*p.E1 + 3/6*p.E6;
    p.E5 = 2/6*p.E1 + 4/6*p.E6;    

    
end
   
    
function     Nd = assign_nodal_coordinates(p, MP)
        
    nBM = MP.NX + 1;

    if isfield(MP,'nz')
        nz = MP.nz;
    else
        nz = round(MP.length_BM/MP.dZ) + 1;
    end
    
    Nd.X = zeros(nz,nBM+13);
    Nd.Y = zeros(nz,nBM+13);
    Nd.Z = zeros(nz,nBM+13);
    
    nBMp = round((1-MP.DC_root_location)*MP.NX)+1;
    nBMa = round(MP.OPC_root_location*MP.NX)+1;
    nBMo = MP.NX-(nBMp-1)-(nBMa-1)+1;
        
    for iz=1:nz    
        nodes_arcuate = linspace(p.A0(iz,1),p.A1(iz,1),nBMa); nodes_arcuate(end)=[];
        nodes_IPC_DC_roots = linspace(p.A1(iz,1),p.A2(iz,1),nBMo);
        nodes_pectinate = linspace(p.A2(iz,1),p.A3(iz,1),nBMp); nodes_pectinate(1)=[];
        
        Nd.X(iz,:) = [nodes_arcuate, nodes_IPC_DC_roots, nodes_pectinate,...
                    0.5*(p.A0(iz,1)+p.CC(iz,1)),0.5*(p.A1(iz,1)+p.CC(iz,1)), p.DD(iz,1), ...
                        p.CC(iz,1), p.BB(iz,1), p.FF(iz,1), p.E1(iz,1), p.E2(iz,1), p.E3(iz,1), p.E4(iz,1), p.E5(iz,1), p.E6(iz,1),p.PP(iz,1)];
        Nd.Y(iz,:) = [zeros(1,nBM),...
                    0.5*(p.A0(iz,2)+p.CC(iz,2)),0.5*(p.A1(iz,2)+p.CC(iz,2)), p.DD(iz,2), ...
                        p.CC(iz,2), p.BB(iz,2), p.FF(iz,2), p.E1(iz,2), p.E2(iz,2), p.E3(iz,2), p.E4(iz,2), p.E5(iz,2), p.E6(iz,2),p.PP(iz,2)];
    end
   
    z_tilt = -0.25*MP.tilt_OHC*10;      % O--+------X
    
% %     %.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.%    
% %     z_tilt = 0;     % To test no-OHC-tilt condition
% %     %.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.%
    
              
    if MP.dim == 2
        z_tilt = 0;
    end
        
    
    Nd.Z = [zeros(1, nBM), ...
                0, 0, 0, ...
                    0,  z_tilt, z_tilt, z_tilt, z_tilt, z_tilt, z_tilt, z_tilt, z_tilt, z_tilt]; % Tilting of OHC and DC process in longitudinal direction
                
    if isfield(MP,'zz')
        zz = MP.zz;
    else
        zz = (-0.5*MP.length_BM:MP.dZ:0.5*MP.length_BM)';
    end
    
    Nd.Z = ones(size(zz))*Nd.Z;

    for ic = 1:length(zz)
        Nd.Z(ic,:)=Nd.Z(ic,:) + zz(ic);
    end
    
    
    Nd.name = cell(size(Nd.X));

    for ii = 1:size(Nd.X,2)
        if  ii< nBM
            Nd.name(:,ii) = {['A' num2str(ii-1)]};
        elseif ii == nBM
            Nd.name(:,ii) = {'AX'};
        elseif ii == nBM + 1
            Nd.name(:,ii) = {'C1'};
        elseif ii == nBM + 2
            Nd.name(:,ii) = {'C2'};
        elseif ii == nBM + 3
            Nd.name(:,ii) = {'DD'};
        elseif ii == nBM + 4
            Nd.name(:,ii) = {'CC'};
        elseif ii == nBM + 5
            Nd.name(:,ii) = {'BB'};
        elseif ii == nBM + 6
            Nd.name(:,ii) = {'FF'};
        elseif ii == nBM + 7
            Nd.name(:,ii) = {'E1'};
        elseif ii == nBM + 8
            Nd.name(:,ii) = {'E2'};
        elseif ii == nBM + 9
            Nd.name(:,ii) = {'E3'};
        elseif ii == nBM + 10
            Nd.name(:,ii) = {'E4'};
        elseif ii == nBM + 11            
            Nd.name(:,ii) = {'E5'};
        elseif ii == nBM + 12            
            Nd.name(:,ii) = {'E0'};
        elseif ii == nBM + 13            
            Nd.name(:,ii) = {'PP'};
        else
            Nd.name(:,ii) = {'XX'};
        end        
    end
    
    for iz=1:size(Nd.X,1)
        for ix=1:size(Nd.X,2)
            ax = Nd.X(iz,ix)/Nd.X(iz,nBM);
            if ix<nBM && abs(ax-MP.OPC_root_location)<1e-2
                Nd.name(iz,ix) = {'AP'};
            elseif ix<nBM && abs(ax-MP.DC_root_location)<1e-2
                Nd.name(iz,ix) = {'AM'};
            end
        end
    end
    
end
    
% % function Nd = adjust_TC(Nd)
% %     % remove mid-nodes in the IPC and OPC
% %     rnodes = strcmp(Nd.name(1,:),'C1') | strcmp(Nd.name(1,:),'C2');
% %     Nd.X(:,rnodes) = [];
% %     Nd.Y(:,rnodes) = [];
% %     Nd.Z(:,rnodes) = [];
% %     Nd.name(:,rnodes) = [];
% % end
    
    
function Nd = initialize_nodes(Nd)
        
    Nd.N = length(Nd.X);
    Nd.Nc = size(Nd.X,2);
    Nd.Nr = size(Nd.X,1);
        
    Nd.dx = zeros(size(Nd.X));
    Nd.dy = zeros(size(Nd.X));
    Nd.dz = zeros(size(Nd.X));
    Nd.rx = zeros(size(Nd.X));
    Nd.ry = zeros(size(Nd.X));
    Nd.rz = zeros(size(Nd.X));
    
    NDOF = 6;
    Nd.BC = ones(Nd.Nr* Nd.Nc, NDOF);
    Nd.N = Nd.Nr*Nd.Nc;
end