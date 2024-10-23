function El = model_moduli5_2021(El, Nd, MP)
% # 
% # function El = model_moduli5_2021(El, Nd, MP)
% #
% # This was updated from 'model_modulid' to reflect the geometry changes
% # matching OCT anatomy.
% # The main (only) difference is the YM of eBMz.
% #
% #
% #
% # Define cross sections and elastic moduli
% #
% # First written: 2009.06.xx
% # Last modified: 2021.04.18
% #
% # Jong-Hoon Nam at University of Rochester
% # 

    global coef

    eBMa = strcmp('BMa',El.name);   % radial BM fiber layer arcuate zone
    eBMp = strcmp('BMp',El.name);   % radial BM fiber layer BM pectinate zone
    eBMz = strcmp('BMza',El.name) | strcmp('BMzp',El.name);   %   
    
    eOHC = strcmp('OHC',El.name);
    ePC  = strcmp('IPC',El.name) | strcmp('OPC',El.name);
    eIPC  = strcmp('IPC',El.name);
    eOPC  = strcmp('OPC',El.name);
    eDCb = strcmp('DCb',El.name);
    eDCp = strcmp('DCp',El.name);
    eDCz = strcmp('DCz',El.name);
    eDCr = strcmp('DCr',El.name);
    eRLx1 = strcmp('RLx1',El.name);
    eRLz = strcmp('RLz',El.name);
    eRLp = strcmp('RLp',El.name);
    eRLpz = strcmp('RLpz',El.name);
    eTMx1 = strcmp('TMx1',El.name); % TM attachment
    eTMx2 = strcmp('TMx2',El.name); % TM body
    eTMz0 = strcmp('TMz0',El.name);
    eTMz1 = strcmp('TMz1',El.name);
    eTMz2 = strcmp('TMz2',El.name);
    eOHB = strcmp('OHB',El.name); % OHC bundle and linear spring Nd1:bb, Nd2:e1
    exHB = strcmp('xHB',El.name); % Pseudo rod at OHC bundle location
    eANK = strcmp('ANK',El.name); % Pseudo link to conveniantly calculate translational displ.
    eTCz = strcmp('TCz',El.name);
    
    El.type = ones(El.N,1)*0.1; % Default Timoshenko beam: 0; Euler beam: 0.1
    El.type(eRLz) = 0.1;
    El.type(eDCz) = 2;
    El.type(eANK) = 2;        
    El.type(eTMx2) = 0.1;
    El.type( eBMz | eTMz1 | eTMz2 | eTCz | eRLz | eRLp) = 0.1;

    El.hinge = zeros(El.N,2); % initiate hinge condition
    El.hinge(eOHB,2) = 1;

    El.oL = zeros(El.N,1);
    El.L = zeros(El.N,1);
    El.odir = zeros(El.N,3);
    El.dir = zeros(El.N,3);
    
    El = get_vector(El, Nd, 0);
    
    El.dim = zeros(El.N,2);  % dim = [bx, bz], note that x is the primary axis    
    El = set_dimension_gradient(El, Nd, MP);    
    
    % Elements with rectangular cross-section
    eRectangular = eBMa | eBMp | eBMz | eRLx1 | eRLz | eRLp | eRLpz | eTMx1 | eTMx2 | eTMz0 | eTMz1 | eTMz2 | eOHB | exHB | eANK;
    % Elements with circular cross-section
    eCircular = eOHC | ePC | eOPC | eDCb | eDCp | eDCz | eDCr | eTCz;
    
    El.A = zeros(El.N,1);
    El.Iz = zeros(El.N,1);
    El.Iy = zeros(El.N,1);
    
    El.A(eRectangular) = El.dim(eRectangular,1).*El.dim(eRectangular,2);
    El.A(eCircular) = (pi/4)*El.dim(eCircular,1).^2;
    
    El.Iz(eRectangular) = (1/12)*El.dim(eRectangular,1).*El.dim(eRectangular,2).^3;
    El.Iz(eCircular) = (pi/64)*El.dim(eCircular,1).^4;
    
    El.Iy(eRectangular) = (1/12)*El.dim(eRectangular,2).*El.dim(eRectangular,1).^3;
    El.Iy(eCircular) = (pi/64)*El.dim(eCircular,1).^4;
        
    % Assigned these values for every element for initialization, if any
    % element has these values, it means that it was not assigned a proper MP.
    El.YM = 123e6*ones(El.N,1);  
    El.SM = 12.3e6*ones(El.N,1);
    El.nu = 0.3*ones(El.N,1);
    
    El.YM(eBMp) = 2.0e9; El.SM(eBMp) = 1/3*El.YM(eBMp);
    El.YM(eBMa) = 2.0e9; El.SM(eBMa) = 1/3*El.YM(eBMa);
 
            
    El = set_material_gradient(El, Nd, MP); % YM for continuous membrane such as BM, TM, RL          
    
    El.YM(eIPC) = coef.e_IPC*400e6;  El.SM(eIPC) = 0.3*El.YM(eIPC);    
    El.YM(eOPC) = coef.e_OPC*400e6;  El.SM(ePC) = 0.3*El.YM(ePC);
    
       
    El.YM(eOHC) = coef.e_OHC*3*15e3; El.SM(eOHC) = 0.3*El.YM(eOHC);

       
    El.YM(eDCb) = coef.e_DCb*20e3; El.SM(eDCb) = 0.3*El.YM(eDCb);
%     El.YM(eDCb) = coef.e_DCb*800e3; El.SM(eDCb) = 0.3*El.YM(eDCb);
    %El.YM(eDCp) = 3e6; El.SM(eDCp) = 0.3*El.YM(eDCp);     
    %El.YM(eDCp) = 12e6; El.SM(eDCp) = 0.3*El.YM(eDCp);
    El.YM(eDCp) = coef.e_DCp*120e6; El.SM(eDCp) = 0.3*El.YM(eDCp);

    El.YM(eDCz) = coef.e_DCz*0.1e3;
    if El.type(eDCz) == 0
        El.SM(eDCz) = 0.001*El.YM(eDCz);
    else      
        El.SM(eDCz) = 0.3*El.YM(eDCz);
    end
    
        
    El.YM(eDCr) = 0e3;    El.SM(eDCr) = 0.3*El.YM(eDCr);

    % note that three DCs are summed so this YM should be divided by 3 to
    % obtain actual value of a single cell
       
    El.YM(eOHB) = 10e6;
    El.SM(eOHB) = 0.3*El.YM(eOHB);
    
    El.YM(exHB) = 10e9;     % This is independently defined in set_material_gradient(), which is not good.
    El.SM(exHB) = 0.3*El.YM(exHB);
    
    El.YM(eANK) = 10e6;
    El.SM(eANK) = 0.3*El.YM(eANK);
    
    El.rho = 1.0e-3*ones(El.N,1); % nanogram/um^3

    El = lumped_mass(Nd, El,  MP);

%     El.rho(eTMx1 | eTMx2 | eTMz0 | eTMz1 | eTMz2) = 0.5e-3;
    El.rho(eTMx1 | eTMx2 | eTMz1 | eTMz2) = 1e-6;
    El.rho(eTMz0) = 1e-6;

    NDOF = 6;
    El.fi = zeros(El.N,2*NDOF);
    
end % of function model_moduli3()
       
function El = set_dimension_gradient(El, Nd, MP)
    
    qMP = consider_geometrical_lgradient(MP);
        
    eBMa = strcmp('BMa',El.name);   % BM arcuate zone
    eBMp = strcmp('BMp',El.name);   % BM pectinate zone
    eBMz = strcmp('BMza',El.name) | strcmp('BMzp',El.name);
    
    eOHC = strcmp('OHC',El.name);
    ePC  = strcmp('IPC',El.name) | strcmp('OPC',El.name);
    eDCb = strcmp('DCb',El.name);
    eDCp = strcmp('DCp',El.name);
    eDCz = strcmp('DCz',El.name);
    eDCr = strcmp('DCr',El.name);
    eRLx1 = strcmp('RLx1',El.name);
    eRLp = strcmp('RLp',El.name);
    eRLz = strcmp('RLz',El.name);
    eRLpz = strcmp('RLpz',El.name);
    eTMx1 = strcmp('TMx1',El.name);
    eTMx2 = strcmp('TMx2',El.name);
    eTMz0 = strcmp('TMz0',El.name);
    eTMz1 = strcmp('TMz1',El.name);  
%     eTMz2 = strcmp('TMz2',El.name);
    eOHB = strcmp('OHB',El.name);
    exHB = strcmp('xHB',El.name);
    eANK = strcmp('ANK',El.name);
    eTCz = strcmp('TCz',El.name);
    
        
    eBMza = strcmp('BMza',El.name);     % [YJ updated,11/01/13]
    eBMzp = strcmp('BMzp',El.name);     % [YJ updated,11/01/13]
    
    El.dim = zeros(El.N,2);  % dim = [by, bz], note that x is the primary axis
                             %      bz is the thickness in primary bending
                             %      direction
    
    eX = eBMa; bz = 'thick_BMa';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);

    eX = eBMp;  bz = 'thick_BMp';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);    

    eX = eBMzp; bz = 'thick_BMp';
    if mean(El.type(eBMz))>2.5 % if the BMz is a plate
        zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
        zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
        El.dim(eX,1) = MP.dZ;
        El.dim(eX,2) = interp1(zz,qMP.width_BM,zi)/MP.NX;        
        El.thick(eX) = interp1(zz,qMP.(bz),zi); % depth
    else              % elseif the BMz is a beam        
        zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
        zz = linspace(min(zi), max(zi), length(qMP.(bz)))';        
        El.dim(eX,1) = interp1(zz,qMP.width_BM,zi)/MP.NX;  
        El.dim(eX,2) = interp1(zz,qMP.(bz),zi); % depth
    end
    
    eX = eBMza; bz = 'thick_BMa';
    if mean(El.type(eBMz))>2.5 % if the BMz is a plate
        zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
        zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
        El.dim(eX,1) = MP.dZ;
        El.dim(eX,2) = interp1(zz,qMP.width_BM,zi)/MP.NX;        
        El.thick(eX) = interp1(zz,qMP.(bz),zi); % depth
    else              % elseif the BMz is a beam        
        zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
        zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
        El.dim(eX,1) = interp1(zz,qMP.width_BM,zi)/MP.NX;  
        El.dim(eX,2) = interp1(zz,qMP.(bz),zi); % depth
    end    
    
%     mthick_BMz = thickness_BMz(Nd, El, eX, MP);
    
    eX = eOHC;  bz = 'diam_OHC';                          
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);    
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';   
    El.dim(eX,1) = interp1(zz,qMP.(bz),zi);               
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);               
    
    eX = ePC; bz = 'thick_PC';                            
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);    
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';   
    El.dim(eX,1) = interp1(zz,qMP.(bz),zi);               
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);               
    
    
    eX = eDCb;  bz = 'diam_DCb';                          
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);    
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';   
    El.dim(eX,1) = interp1(zz,qMP.(bz),zi);               
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);               
    
    eX = eDCp;  bz = 'diam_DCp';                          
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);    
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';   
    El.dim(eX,1) = interp1(zz,qMP.(bz),zi);                
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);               
    
    eX = eDCz;
    El.dim(eX,1) = MP.diam_DCz;
    El.dim(eX,2) = MP.diam_DCz;
    
    eX = eDCr;
    El.dim(eX,1) = MP.diam_DCz;
    El.dim(eX,2) = MP.diam_DCz;
    
    eX = eRLx1; bz = 'thick_RL';
    zi = Nd.Z(El.Nd2(eX));                              % Nd1 to Nd2, MPN
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);
    
    eX = eRLp; bz = 'thick_RL';
    zi = Nd.Z(El.Nd2(eX));                              % Nd1 to Nd2, MPN
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);
    
    eX = eRLpz; bz = 'thick_RL';
    zi = Nd.Z(El.Nd2(eX));                              % Nd1 to Nd2, MPN
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);

    eX = eRLz; bz = 'thick_RL';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);    
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
%     El.dim(eX,1) = MP.dX;                             % no ground for this value
    El.dim(eX,1) = interp1(zz,qMP.width_BM,zi)/MP.NX;  
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);
    
    eX = eTMx1; bz = 'thick_TMa';   % TM attachement thickness
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = MP.dZ;    
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);
    % attachement region has half the thickness of the body
    
    eX = eTMx2; bz = 'thick_TMb';   % TM body thickness
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = MP.dZ;
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);

    eX = eTMz0; bz = 'thick_TMb';                       % Longitudianl TM above the HB
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = 1/6*interp1(zz,qMP.width_TM,zi);     % as if there be 1/12 overhang,[YJ updated,10/24/13]
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);

    eX = eTMz1; bz = 'thick_TMb';                       % Longitudianl TM at body
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = 1/6*interp1(zz,qMP.width_TM,zi);    
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);
    
    
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
    eX = eOHB;    
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    bz = 'k_OHB';   kHB = interp1(zz,qMP.(bz),zi);
    bz = 'height_HB';   hHB = interp1(zz,qMP.(bz),zi);
    YM_HB = 10e6;
    a1 = 0.9;
        
    wHB = 1.0;  % HB width [um]
    El.dim(eX,1) = wHB;
    El.dim(eX,2) = hHB.*( (4*a1*kHB/YM_HB/wHB).^(1/3));        % Cantilever beam stiffness kHB = 3*YM*Iz/H^3, Iz = (w*t^3)/12 --> t = H*(4*kHB/YM/w)^1/3
    
    eX = exHB;
    El.dim(eX,1) = wHB;
    El.dim(eX,2) = 0.5;
    
    eX = eANK;
    El.dim(eX,1) = wHB;
    El.dim(eX,2) = (1-a1)*kHB.*El.oL(eX)/YM_HB/wHB;            % Bar stiffness kANK = YM*w*t/L, where L = El.oL --> t = kANK*L/YM/w
    
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%
%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%%N%    
    
    eX = eTCz; bz = 'thick_TCz';
    zi = mean([Nd.Z(El.Nd1(eX)),Nd.Z(El.Nd2(eX))], 2);
    zz = linspace(min(zi), max(zi), length(qMP.(bz)))';
    El.dim(eX,1) = interp1(zz,qMP.(bz),zi);
    El.dim(eX,2) = interp1(zz,qMP.(bz),zi);
    
end % of function set_dimension_gradient()
    
function  El = set_material_gradient(El, Nd, MP) 
% YM for continuous membrane such as BM, TM, RL

    qMP = consider_MP_lgradient(MP);    
    eBMz = strcmp('BMza',El.name) | strcmp('BMzp',El.name);
    eRLx1 = strcmp('RLx1',El.name);
    eRLp = strcmp('RLp',El.name);
    eRLz = strcmp('RLz',El.name);
    eTMx1 = strcmp('TMx1',El.name); 
    eTMx2 = strcmp('TMx2',El.name);
    eTMz = strcmp('TMz0',El.name) | strcmp('TMz1',El.name);
    eTCz = strcmp('TCz',El.name);
    
    zi = mean([Nd.Z(El.Nd1),Nd.Z(El.Nd2)], 2);
    
    mp_YM = {'YM_BMz'; 'YM_RLx1'; 'YM_RLp'; 'YM_RLz'; 'YM_TCz'; 'YM_TMx1'; 'YM_TMx2'; 'YM_TMz'};
    mp_SM = {'SM_BMz'; 'SM_RLx1'; 'SM_RLp'; 'SM_RLz'; 'SM_TCz'; 'SM_TMx1'; 'SM_TMx2'; 'SM_TMz'};
    lg_elems = [eBMz, eRLx1, eRLp, eRLz, eTCz, eTMx1, eTMx2, eTMz];

    for kk=1:size(lg_elems,2)
        pYM = mp_YM{kk}; 
        pSM = mp_SM{kk}; 
        elem = lg_elems(:,kk);
        
        zk = linspace(min(zi(elem)), max(zi(elem)), length(qMP.(pYM)))';
            
        El.YM(elem) = interp1(zk, qMP.(pYM), zi(elem));
        El.SM(elem) = interp1(zk, qMP.(pSM), zi(elem));
    end
    
end % of function set_material_gradient



function C = consider_MP_lgradient(MP)

        
    A = MP_of_longitudinal_elements('a');
    B = MP_of_longitudinal_elements('b');

    zz = (-0.5*MP.length_BM:MP.dZ:0.5*MP.length_BM)';
    if abs(MP.opt_lgradient) >= 1       % no gradient, the value at z = MP.opt_lgraident is used for the entire span
        zz = MP.opt_lgradient*zz;       % For example,  lgradient=1, forward gradient; lgradient=-1 means reverse gradient
    else
        zz = median(zz)*ones(size(zz));
    end
    

    if MP.loc<=6       % values will be referenced to the basal properties at Z = 2mm
        zz = 1e3*(MP.loc-2) + zz;
        C0 = B;
    else                % values will be referenced to the apiaal properties at Z = 10mm
        zz = 1e3*(MP.loc-10) + zz;
        C0 = A;
    end
    
    
    POC = fieldnames(A);
       
    for ii=1:length(POC)
        prop = POC{ii};
        
        %grad = set_longi_gradient(A,B,prop,pscal);
        %C.(prop) = MP.(prop)*exp(grad*zz);
        grad = set_longi_gradient(A,B,prop,MP.gscale);

        if strcmpi(MP.gscale(1:3),'LIN')
            C.(prop) = C0.(prop) + grad*zz;
        elseif strcmpi(MP.gscale(1:3),'LOG')
            C.(prop) = C0.(prop)*exp(grad*zz);
        else
            C.(prop) = C0.(prop) + 0*zz;
        end

        if min(C.(prop))<0
            wstring = {'MP gradient inappropriate';...
                ['Resuring in negative value of [' prop ']']; ...
                'Gradient ignored'};
            wname = 'At consider_MP_lgradient.m';
            warndlg(wstring, wname);
            C.(prop) = mean(C.(prop))*ones(size(zz));
        end        
    end        
    
end



function C = MP_of_longitudinal_elements(loc)

    global coef

    %C.YM_BMz = 400e3;    C.SM_BMz = 1/3*C.YM_BMz;
    C.YM_BMz = coef.e_BMz*100e3;    C.SM_BMz = 1/3*C.YM_BMz;
%     C.YM_BMz = coef.e_BMz*100e3;    C.SM_BMz = 0.001*C.YM_BMz;

    C.YM_TCz = coef.e_YMz*0.2e6;    C.SM_TCz = 0.3*C.YM_TCz;
%     C.YM_TCz = coef.e_YMz*0.2e6;    C.SM_TCz = 0.001*C.YM_TCz;

    C.YM_TMz = coef.e_TMz*2e3;   C.SM_TMz = 0.3*C.YM_TMz;   % A factor of 10 to match Freeman et al.'s TM traveling wave form
%     C.YM_TMz = coef.e_TMz*2e3;   C.SM_TMz = 0.001*C.YM_TMz;
    if strcmpi(loc,'b')
    % Material properties at 2 mm from BASE to apex
    %%%%%% Low longitudinal couple        
     
%         C.YM_RLx1 = coef.e_RLx1*5*100e6;   C.SM_RLx1 = 0.3*C.YM_RLx1;        
        C.YM_RLx1 = coef.e_RLx1*50e6;   C.SM_RLx1 = 0.3*C.YM_RLx1; 

        C.YM_RLp = coef.e_RLp*20e6;   C.SM_RLp = 0.3*C.YM_RLp;
%         C.YM_RLp = coef.e_RLp*20e6;   C.SM_RLp = 0.001*C.YM_RLp;

        C.YM_RLz = coef.e_YMz*0.2e6;    C.SM_RLz = 0.3*C.YM_RLz;
%         C.YM_RLz = coef.e_YMz*0.2e6;    C.SM_RLz = 0.001*C.YM_RLz;

        C.YM_RLpz = coef.e_YMz*0.2e6;    C.SM_RLz = 0.3*C.YM_RLpz;
        
        C.YM_TMx1 = coef.e_TMx1_B*100e3;   C.SM_TMx1 = 0.3*C.YM_TMx1;     % stiff values
        C.YM_TMx2 = 1/4*C.YM_TMx1; C.SM_TMx2 = 0.3*C.YM_TMx2;

    elseif  strcmpi(loc,'a')
        % Material properties at 10  mm from base to APEX
        
%         C.YM_RLx1 = coef.e_RLx1*5*100e6;   C.SM_RLx1 = 0.3*C.YM_RLx1;        
        C.YM_RLx1 = coef.e_RLx1*50e6;   C.SM_RLx1 = 0.3*C.YM_RLx1; 
        C.YM_RLp = coef.e_RLp*5.0e6;   C.SM_RLp = 0.3*C.YM_RLp;
%         C.YM_RLp = coef.e_RLp*5.0e6;   C.SM_RLp = 0.001*C.YM_RLp;

        C.YM_RLz = coef.e_YMz*0.2e6;    C.SM_RLz = 0.3*C.YM_RLz;
%         C.YM_RLz = coef.e_YMz*0.2e6;    C.SM_RLz = 0.001*C.YM_RLz;

        C.YM_RLpz = coef.e_YMz*0.2e6;    C.SM_RLz = 0.3*C.YM_RLz;       
                
        C.YM_TMx1 = coef.e_TMx1_A*20e3;   C.SM_TMx1 = 0.3*C.YM_TMx1;
%         C.YM_TMx1 = coef.e_TMx1_A*100e3;   C.SM_TMx1 = 0.3*C.YM_TMx1;
        C.YM_TMx2 = 1/4*C.YM_TMx1; C.SM_TMx2 = 0.3*C.YM_TMx2;       
        
    else
    
    end   
end
    