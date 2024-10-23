function [C,G,V0] = mElectrical2_cmod_051216(MP,OHC_sg)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This is the analysis of whole three scalae circuit.This is giving output for dynamic condition.
% % This circuit is for solving the differential equation for the circuit formed in the function "tj_OClongnew_24thapril_real_dynamic2"
% % No input for this function explicitly but it reads the BM vibration from the mat file 'input', the output is time and node voltages
% % First written in April,2013 as a part of the course project BME 404
% % The circuit is a combination Mistrik et al. 2009 paper and NAMLAB circuit

   % Initial condition
   
   nz = round(MP.length_BM/MP.dZ)+1;           % # of sections along longitudinal direction
   nv = 5;
   % MSA: these definition for nodes are not clear in where they sit in the
   % matrix for the equations
    % % Voltage vector, V= [ V1,1; V1,2; V1,4; V1,6; V1,7; V1,8; V1,9];
    % % V1,1= Node of scala media
    % % V1,2= Intracellular space of OHC
    % % V1,4= Extracellular space of OHC
    % % V1,6= Intracellular space of IHC
    % % V1,7= Extracellular space of IHC
    % % V1,8= Node of scala tympani
    % % V1,9= Node of scala vestibuli
    % % %%%%%%%%Read input file
    [StV,OHC,IHC,SL,OoC] = organ_Corti_parameters_wgrad(nz,OHC_sg);

    [G0,f]= make_G0matrix(nv,nz,StV,OHC,IHC,SL,OoC);
    V0 = G0\f;
    [C,G] = make_CGmatrix(nv,nz,V0,StV,OHC,IHC,SL,OoC,OHC_sg);

end

function [C,G] = make_CGmatrix(nv,nz,V0,StV,OHC,IHC,SL,OoC,OHC_sg)

    C=zeros(nv*nz,nv*nz);


    for iz=1:nz
        M = OHC_sg(iz).M;
        Vm0 = V0((iz-1)*5+2)-V0((iz-1)*5+3);
        [a0,Cm0] = coeff_Vm(M,Vm0);
        OHC.Cm(iz) = Cm0;

        Cm=[OHC.Cs(iz)+StV.C(iz)+IHC.Cs(iz) -OHC.Cs(iz) 0 -IHC.Cs(iz) 0;
            -OHC.Cs(iz) OHC.Cs(iz)+OHC.Cm(iz)+a0 -OHC.Cm(iz)-a0 0 0;
            0 -OHC.Cm(iz)-a0 OHC.Cm(iz)+a0 0 0;
            -IHC.Cs(iz) 0 0 IHC.Cs(iz)+IHC.Cm(iz) -IHC.Cm(iz);
            0 0 0 -IHC.Cm(iz) IHC.Cm(iz)];
        C(5*iz-4:5*iz,5*iz-4:5*iz)=Cm;
    end
    C = sparse(C);


    % Forming the matrices
    % [StV,OHC,IHC,SL,OoC]=organ_Corti_parameters_wgrad(n,yBM);
    % [SV,SM,ST]=scala_parameters(n);

    G=zeros(nv*nz,nv*nz);
    % Connectivity matrix
    for iz=1:nz
        i=5*iz-4;
        if iz~=1
         G(i,i-5)=-StV.G2(iz);
         G(i+2,i-3)=-OoC.G2(iz);
         G(i+4,i-1)=-SL.G2(iz);
        end
        if iz~=nz
          G(i,i+5)=-StV.G2(iz);
         G(i+2,i+7)=-OoC.G2(iz);
         G(i+4,i+9)=-SL.G2(iz);
        end
    end
    
    for iz=1:nz
        Vm0 = V0((iz-1)*nv+2)-V0((iz-1)*nv+3);
        M = OHC_sg(iz).M;
        [Gm0,sGm] = coeff_OM(M,Vm0);
        az = Gm0+sGm*75; % MSA: hard-coded value!!!!
        G1=[StV.G1(iz)+StV.G2(iz)+OHC.Gs(iz)+IHC.Gs(iz) -OHC.Gs(iz) 0 -IHC.Gs(iz) 0 ];
        G2=[-OHC.Gs(iz) OHC.Gs(iz)+az -az 0 0];
        G3=[0 -az OoC.G1(iz)+az+OoC.G2(iz) 0 0];
        G4=[-IHC.Gs(iz) 0 0 IHC.Gs(iz)+IHC.Gm(iz) -IHC.Gm(iz)];
        G5=[0 0 0 -IHC.Gm(iz) SL.G1(iz)+SL.G2(iz)+IHC.Gm(iz)];
        Gm=[G1;G2;G3;G4;G5];
        if iz==1 || iz==nz
            G((5*iz-4):5*iz,(5*iz-4):5*iz)=Gm;
        else
            G((5*iz-4):5*iz,(5*iz-4):5*iz)=Gm+diag([StV.G2(iz) 0 OoC.G2(iz) 0 SL.G2(iz)]);
        end
    end
    G = sparse(G);

end

function [Gm0,sGm] = coeff_OM(M,Vm0)

    zeta0 = (1+exp(-(Vm0-M.VG05)/M.eK))^(-1);
    Gm0 = M.Gmax*zeta0;
    sGm = -M.Gmax*zeta0^2*exp(-(Vm0-M.VG05)/M.eK)*(-1/M.eK);

end

function [StV,OHC,IHC,SL,OoC] = organ_Corti_parameters_wgrad(nz,OHC_sg)
    % d=0.1;    %source: Mistrik 2009
    % BM.x=yBM;    % Basilar membrane displacement (input)
    
    % Stria Vascularis(StV)
    StV.C = 1/4*(66.7)*ones(nz,1);  % Stria Vascularis capacitance in PF, chosen from Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    
    StV.G1 = 1/4*(6.67*1e3*5)*ones(nz,1);     
    StV.G2 = (1.7*1e4)*ones(nz,1); %%%%changed on 08/14/19
        
    StV.E= 90*ones(nz,1);
        
    % Outer Hair Cell(OHC)
    for iz=1:nz
        % source: Yanju's single cell capacitance
        OHC.Cs(iz) = 3*OHC_sg(iz).S.C; % OHC steriocilial capacitance in pF,   All OHC parameters is referred from Yanju and Nam 2012
        OHC.Cm(iz) = 3*OHC_sg(iz).M.CLin; % OHC basolateral linear capacitance in pF
        OHC.Gs_max(iz) = 3*OHC_sg(iz).S.Gmax;%OHC steriocilial conductance(maximum)in nS, %MSA: temporarily changed conductance by a factor to later figure out why it is so high in the model
        OHC.Gmax(iz) = 3*OHC_sg(iz).M.Gmax;%OHC basolateral conductance in nS %MSA: temporarily changed conductance by a factor to immitate experiment
        OHC.Ek(iz) = 75; %OHC equilibrium potential in mV, source: Mistrik 2009
        %%% incorporating open probability
        OHC.Gs(iz)= OHC.Gs_max(iz).*OHC_sg(iz).ini.po;
        OHC.ini.zeta(iz) = OHC_sg(iz).ini.zeta;
    end
    
    % Inner Hair Cell (IHC)
    IHC.Cs = ones(nz,1);% OHC steriocilial capacitance in pF, source: Mistrik 2009
    IHC.Cm = 10*ones(nz,1); % OHC basolater linear capacitance in pF, source: Mistrik 2009
    IHC.Gs_max = 0.25*(99.5)*ones(nz,1); %IHC steriocilial conductance(maximum)in nS, source: Mistrik 2009, "The deviant factor is 0.6 "
    IHC.Gm = 0.25*(33.333)*ones(nz,1);        %IHC basolateral conductance in nS, source: Mistrik 2009
    IHC.Ek = 75*ones(nz,1); %IHC equilibrium potential in mV, source: Mistrik 2009
    %     IHC.Cs = 0*ones(n,1);% OHC steriocilial capacitance in pF, source: Mistrik 2009
    %     IHC.Cm = 0*ones(n,1); % OHC basolater linear capacitance in pF, source: Mistrik 2009
    %     IHC.Gs_max = 0*(99.5)*ones(n,1); %IHC steriocilial conductance(maximum)in nS, source: Mistrik 2009
    %     IHC.Gm = 0*(33.333)*ones(n,1);        %IHC basolateral conductance in nS, source: Mistrik 2009
    %     IHC.Ek = 0*ones(n,1); %IHC equilibrium potential in mV, source: Mistrik 2009
    % poi = 0.02;
    poi=0.1;
    IHC.Gs = IHC.Gs_max*poi;
    
    % Spiral Limbus (SL)
    SL.G1 = 0.25*(99.5)*ones(nz,1); % Spiral limbus conductance in nS, source : Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    SL.G2 = 0.25*(33.33)*ones(nz,1);     % Spiral limbus connectivity conductance in nS, source : Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    %     SL.G1 = 0*(99.5)*ones(n,1); % Spiral limbus conductance in nS, source : Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    %     SL.G2 = 0*(33.33*45)*ones(n,1);     % Spiral limbus connectivity conductance in nS, source : Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    % organ of Corti (OoC)
    OoC.G1 = 0.25*6.67*1e4*ones(nz,1);   % Organ_of_corti conductance in nS, source: Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5) %MSA: temporarily changed conductance by a factor to later figure out why it is so high in the model
    OoC.G2 = 0.25*6.67*1e4*ones(nz,1); % Organ_of_corti connectivity conductance in nS, source: same as OoC.G1 for Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    %     OoC.G1 = 10000000*(6.67e3)*ones(n,1);   % Organ_of_corti conductance in nS, source: Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
    %     OoC.G2 = 0*(1905*7)*ones(n,1); % Organ_of_corti connectivity conductance in nS, source: Mistrik 2009,modified for 40 micron instead of 60 micron(multiplied by 1.5)
end

function [alpha0,Cm0,M] = coeff_Vm(M,Vm0)

%     a = 37.4e-3; % e/kT (1/mV),Note the unit conversion factor 1e-3, V -> mV
%     epsilon = exp(-M.zC*a*(Vm0-M.VC05));
% 
%     Cm0 = M.Qmax*(M.zC*a)*epsilon/(1+epsilon)^2 + M.CLin;
% 
%     alpha0 = -Vm0*M.Qmax*(M.zC*a)^2*(epsilon/(1+epsilon)^2 - 2*epsilon^2/(1+epsilon)^3);
% ****************************************************************************************
% Unit conversion
% Cm: pF
% Gm: nS
% Q: pC
    a = 37.4; % e/kT (1/V),
    epsilon = exp(-M.zC*a*(Vm0-M.VC05)*1e-3);%Note the unit conversion factor 1e-3, V -> mV

    Cm0 = M.Qmax*(M.zC*a)*epsilon/(1+epsilon)^2 + M.CLin;

    alpha0 = -1e-3*Vm0*M.Qmax*(M.zC*a)^2*(epsilon/(1+epsilon)^2 - 2*epsilon^2/(1+epsilon)^3);%Note the unit conversion factor 1e-3, V -> mV

end

function [G,f] = make_G0matrix(nv, nz,StV,OHC,IHC,SL,OoC)
    % % This is the analysis of whole three scalae circuit.This is giving output for dynamic condition.
    % % This circuit is for creating the differential equation for the circuit
    % % First written in May,2014 as a part of the course project BME 404
    % % The circuit is a combination Mistrik et al. 2009 paper and NAMLAB circuit
    % % The circuit takes time, voltage and basilar membrane displacement as input and returns differential voltage as output

    % Forming the matrices
    % [StV,OHC,IHC,SL,OoC]=organ_Corti_parameters_wgrad(n,yBM);
    % [SV,SM,ST]=scala_parameters(n);

    
    G=zeros(nv*nz,nv*nz);
    % MSA: this entire code needs further explanation. The vector of
    % variables should be defined, the different matrices should be
    % explained and what is coded here is different from thesis
    % Connectivity matrix
    for iz=1:nz
        i=5*iz-4;
        if i~=1
            G(i,i-5)=-StV.G2(iz);
            G(i+2,i-3)=-OoC.G2(iz);
            G(i+4,i-1)=-SL.G2(iz);
        end
        if i~=(5*nz-4)
            G(i,i+5)=-StV.G2(iz);
            G(i+2,i+7)=-OoC.G2(iz);
            G(i+4,i+9)=-SL.G2(iz);
        end  
    end
    for iz=1:nz
        fs=[StV.E(iz)*StV.G1(iz);-OHC.Ek(iz)*OHC.ini.zeta(iz)*OHC.Gmax(iz);OHC.Ek(iz)*OHC.ini.zeta(iz)*OHC.Gmax(iz);-IHC.Ek(iz)*IHC.Gm(iz);IHC.Ek(iz)*IHC.Gm(iz)];
      f(5*iz-4:5*iz,1)=fs;
    end
    for iz=1:nz
        G1=[StV.G1(iz)+StV.G2(iz)+OHC.Gs(iz)+IHC.Gs(iz) -OHC.Gs(iz) 0 -IHC.Gs(iz) 0];
        G2=[-OHC.Gs(iz) OHC.Gs(iz)+OHC.ini.zeta(iz)*OHC.Gmax(iz) -OHC.ini.zeta(iz)*OHC.Gmax(iz) 0 0];
        G3=[0 -OHC.ini.zeta(iz)*OHC.Gmax(iz) OoC.G1(iz)+OHC.ini.zeta(iz)*OHC.Gmax(iz)+OoC.G2(iz) 0 0];
        G4=[-IHC.Gs(iz) 0 0 IHC.Gs(iz)+IHC.Gm(iz) -IHC.Gm(iz)];
        G5=[0 0 0 -IHC.Gm(iz) SL.G1(iz)+SL.G2(iz)+IHC.Gm(iz)];
        Gm=[G1;G2;G3;G4;G5];
        if iz==1 || iz==nz
            G((5*iz-4):5*iz,(5*iz-4):5*iz)=Gm;
        else
            G((5*iz-4):5*iz,(5*iz-4):5*iz)=Gm+diag([StV.G2(iz) 0 OoC.G2(iz) 0 SL.G2(iz)]);
        end
    end
    G = sparse(G);
end
