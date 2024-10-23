function compute_relative_stiffness(MP,Nd,El)

    if ~exist('MP','var')
%         set_Global_coefficients; % this is setting coefficients for parameter study
        
        dim = 3;                    % dimension 2 or 3-D
        opt_lgradient = 1;          % gradient along the longitudinal direction
        opt_MP_set = 1;         % 0: the MP set before 2021.04
    
        MP.loc = 6;             % [mm]
        MP.length_BM = 12e3;   % [um]      
        MP.xx = MP.loc-0.5*MP.length_BM*1e-3:0.01:MP.loc+0.5*MP.length_BM*1e-3;
        
        [Nd, El, MP] = model_3D(MP,dim,opt_lgradient, opt_MP_set);
    end

    [~,kTMa] = compute_TM_axial_stiffness(Nd, El);
    [~,kTMb] = compute_TM_bending_stiffness(Nd, El);
    [kOHB_a,kOHB_b] = compute_stiffness(MP,El,'ANK');
    [kOHC_a,kOHC_b] = compute_stiffness(MP,El,'OHC');
    [kDCb_a,kDCb_b] = compute_stiffness(MP,El,'DCb');
    [kRLx1_a,kRLx1_b] = compute_stiffness(MP,El,'RLx1');

    figure; clf;
    semilogy(MP.xx,kTMa); hold on;
    semilogy(MP.xx,kTMb); hold on;
    semilogy(MP.xx,kOHB_a); hold on;
    semilogy(MP.xx,kOHB_b); hold on;
    semilogy(MP.xx,kOHC_a); hold on;
    semilogy(MP.xx,kOHC_b); hold on;
    semilogy(MP.xx,kDCb_a); hold on;
    semilogy(MP.xx,kDCb_b); hold on;
    semilogy(MP.xx,kRLx1_a); hold on;
    semilogy(MP.xx,kRLx1_b); hold on;    
    xlim([0 12]);

    naxis;
    legend('kTMa','kTMb','kOHB_a','kOHB_b','kOHC_a','kOHC_b','kDCb_a','kDCb_b','kRLx1_a','kRLx1_b');

end

function [ka, kb] = compute_stiffness(MP,El,eName)

    id = find(strcmp(El.name,eName));
    nn = 1:length(id)/length(MP.xx):length(id);

    YM = El.YM(id(nn));
    L = El.oL(id(nn));
    A = El.A(id(nn));
    Iz = El.Iz(id(nn));


    ka = YM.*A./L;
    kb = 3*YM.*Iz./L.^3;

    if strcmp(eName,'ANK')
        ka = 10*ka;
    end

    ka = ka*1e-3; % from [uN/m] to [mN/m]
    kb = kb*1e-3; % from [uN/m] to [mN/m]

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