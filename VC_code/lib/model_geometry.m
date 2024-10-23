function MP = model_geometry(MP, loc, opt_lgradient)
% # function MP = model_geometry(MP, loc)
% #
% # Assign geometrical properties
% #
% # First written: 2009.06.xx
% # Last modified: 2009.11.25
% # Written by: Jong-Hoon Nam at University of Rochester
% #
% #  
    
    
    if strcmpi(loc(1),'A'),   % Apical turn
    % Apical part 'A' targets 10 mm from stapes
       
        MP.k_OHB = 3*3e3; %3*4.0e3; % outer hair cell bundle stiffness, 5 mN/m  times 3 rows of OHC       
      
        MP.thick_BMp = 0.7; % BM thickness at pectinate zone        
        
        
        MP.thick_BMa = 0.25*MP.thick_BMp; % BM thickness at arcuate zone        
        MP.thick_fluid = 0; % added fluid mass, this is an ad hoc parameter for OCC FE analysis only
        MP.thick_BMz = 55; % 55.5 um Edge et al., 1998, Table 2, BM(f)
        
        MP.rOCmass = 0.5; % Additinoal OC mass 
        
        MP.opt_BM = 0.1;    % BMz element type. 0: Euler beam, 0.1: Timoshenko beam, 2: link
        MP.width_BM = 280;  % radial width                        
        
        %MP.length_BM = MP.length_BM; % length, defined in themain funcation 'virtualCochcleaHarmonicXX'.
            
        MP.thick_RL = 1;  % RL thickness        
        MP.thick_TCz = 5;   % thickness of longitudial beam through TC tips        
        
        MP.diam_DCb = 10;    % Deiter's cell diameter
        MP.diam_DCp = 1;    % Phalangeal process diamter
        MP.diam_DCz = 1;   % Deiter's cell longitudinal coupler diameter
  
        
        MP.thick_TMa = 25;   % thickness of TM at attachement  
        
%         MP.thick_TMa = 15;   % thickness of TM at attachement  
        
        
        MP.thick_TMb = 50;   % thickness of TM at body, half this value at attachment
        MP.thick_TMza = MP.thick_TMa;   % thickness of TM at body, half this value at attachment
        MP.thick_TMzb = MP.thick_TMb;   % thickness of TM at body, half this value at attachment
        
        MP.width_TM = 200;  % width TM, Edge et al., 1998, L_BMp = 200 um, L_TM = 210, L_BMa ~ L_OPC = 110
        MP.height_TMSLa = 65; % height of TM-SLa attachment point
               
        MP.thick_PC = 2;    % diameter of pillar cell, Kikuchi et al., 1991 inferred from # of microtubules

        MP.width_TC = 21;
        MP.height_TC = 75;
        MP.height_HB = 6;
        MP.angle_RC = 20;   % angle of RL wrt BM in degrees, Edge et al., 1998
        
        
        MP.diam_OHC = 7.5;  % 7.5 um He et al, Hear. Res. 1994        
        MP.tilt_OHC = 3;    % # of cell spacing between OHC and DC process tops    
        
        
        
        
        MP.diam_OHC = 9.0;  % 7.5 um He et al, Hear. Res. 1994        
              
        

        
        MP.DD_alpha = 0.6;  % The location of OHC-DC joint is alpha from BB towards AM, see model_nodes.m
        MP.DD_theta = 10;   % The radial tilt angle of OHC
         
        MP.OPC_root_location = 1/3;
        MP.DC_root_location = 1/2;

        
        MP.NX = 12;
        MP.dZ = 10;         % OHC longitudinal spacing 
        MP.NZ = round(MP.length_BM/MP.dZ);
                
    elseif strcmpi(loc(1),'B'), % Basal turn
    % Basal part 'B' targets 2 mm from stapes

% %         MP.k_OHB = 3*50.0e3; % 50 pN/nm outer hair cell bundle stiffness;                 
% %         MP.thick_BMp = 3.5;  % BM thickness at pectinate zone, 3.0 um according to Schweizer et al., 1996 
        MP.k_OHB = 3*40e3; % 40 pN/nm outer hair cell bundle stiffness times 3 rows of OHC;
        
            
        MP.thick_BMp = 3.2;  % BM thickness at pectinate zone, 3.0 um according to Schweizer et al., 1996 
        MP.thick_BMa = 0.25*MP.thick_BMp;  % zona acuata is more compliant
% %         MP.thick_BMa = MP.thick_BMp;  % zona acuata is more compliant
        
        MP.thick_fluid = 0; % added mass effect due to fluid, this is an ad hoc parameter for OCC FE analysis only
        MP.thick_BMz = 35; % 35; % 33.4 um Edge et al., 1998, Table 2, BM(f)
                                % a factor of 1.5 was multiplied to
                                % consider the mass of OC.
        
        MP.rOCmass = 0.5; % Additinoal OC mass 
        
        MP.opt_BM = 0.1;    % BMz element type, 0.0 Euler beam, 0.1 Timoshenko beam
        MP.width_BM = 160;  % radial width

        %MP.length_BM = MP.length_BM; % length, defined in themain funcation 'virtualCochcleaHarmonicXX'.
        MP.thick_RL = 2;    % RL thickness        
        MP.thick_TCz = 5;   % thickness of longitudial beam through TC tips
        
        MP.diam_DCb = 10;    % Deiter's cell diameter
        MP.diam_DCp = 1.5;    % Phalangeal process diamter    
        MP.diam_DCz = 1.0;   % Deiter's cell longitudinal coupler diameter
                

% %       % This properties are closer to Edge et al., 1998 report.
% %         MP.thick_TMa = 20;   % thickness of TM at attachement
        
        MP.thick_TMa = 20;   % thickness of TM at attachement        
        MP.thick_TMb = 30;   % thickness of TM at body, half this value at attachment        
        
        MP.thick_TMza = MP.thick_TMa;   % thickness of TM at attachment
        MP.thick_TMzb = MP.thick_TMb;   % thickness of TM at body
        
        MP.width_TM = 100; %100;      % width of TM, Edge et al., 1998

        MP.height_TMSLa = 75;   % height of TM-SLa attachment point
                           
        MP.thick_PC = 8;    % diameter of pillar cell        
                        
        MP.width_TC = 21;
        MP.height_TC = 46;
        
        MP.height_HB = 2.0;        
        MP.angle_RC = 10;   % angle of RL wrt BM in degrees, Edge et al., 1998                
        
        
        
        MP.diam_OHC = 9.0;  % 9 um He et al, Hear. Res. 1994
        
        
        
        
        %MP.diam_OHC = 16;  % test
        
        
                
        
        MP.tilt_OHC = 3;    % tilt_OHC=2, two cell spacing between OHC and DC process tops
        

        MP.DD_alpha = 0.7;  % The location of OHC-DC joint is alpha from BB towards AM, see model_nodes.m
        MP.DD_theta = 5;   % The radial tilt angle of OHC
        
        MP.OPC_root_location = 1/3; % 1/3 BM width from the inner edge
        MP.DC_root_location = 1/2;  % 1/2 BM width from the inner edge
        MP.NX = 12;         % number of division of radial BM elements
        
        
        MP.dZ = 10;         % OHC longitudinal spacing 
        MP.NZ = round(MP.length_BM/MP.dZ);        

    else  
        % assign MP defending on given tonotopic location
    end
    
    MP.opt_lgradient = 1;
    if exist('opt_lgradient','var'),        
        if isnumeric(opt_lgradient),
            MP.opt_lgradient = opt_lgradient;
        end        
    end
    
    MP.gscale = 'log';       % geometry gradient scale either 'log' or 'linear';
   
end