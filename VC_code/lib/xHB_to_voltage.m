function [xx_t,xa_t,po_t,voltage_t,tt] = xHB_to_voltage(loc,freq,amplitude)

    fftw('dwisdom',[]);
    fftw('planner','patient');

    Tfin = 4/freq; % [ms]
    sRate = freq*40; % [kHz]
    dt = 1/sRate; % [ms]    
    tt = 0:dt:Tfin-dt; % [ms]
    nt = length(tt);
    ff = (0:nt-1)/Tfin;
    ffsym = (-nt/2:nt/2-1)/Tfin;

    MP.length_BM = 0;
    MP.dZ = 10;
    MP.loc = loc;
    MP.freq = freq;
    MP.d = amplitude;

    %% model OHC circuit
    nOHC = 1;
    nv = 5;

    OHC(nOHC) = struct('HB',[],'S',[],'C',[]);

    za = 10; zb = 2;        
    if MP.loc<=6
        zz = 1e3*(MP.loc-zb);
    else
        zz = 1e3*(MP.loc-za);
    end

    for iz = 1:nOHC

        [S, M] = OHC_electrical_params2(MP.loc,zz(iz));

        HB = HB_transduction_params2(MP.loc, zz(iz));
        HB = ini_OHC(S,HB);
        
        OHC(iz).count = 3;
        OHC(iz).HB = HB;
        OHC(iz).S = S;
        OHC(iz).M = M;
    
        OHC(iz).S.po0 = HB.po0;
        OHC(iz).ini.po = HB.po0;
        OHC(iz).ini.Gs = OHC(iz).ini.po*OHC(iz).S.Gmax;
        OHC(iz).ini.zeta = 0.8;
        OHC(iz).M.zeta = 0;

    end

    IHC = model_IHC(nOHC);
    err = 1;

%     while err > 1e-12
        [~,~,~,Vnd0] = mElectrical2_cmod_MSA(MP,OHC,IHC);
        M = [OHC.M];
        ini = [OHC.ini];
        Vnd0 = reshape(Vnd0,nv,nOHC);
        Vm0 = Vnd0(2,:) - Vnd0(3,:);
        zeta0 = (1+exp(-(Vm0-[M.VG05])./[M.eK])).^(-1);
%         err = norm(zeta0 - [ini.zeta])/norm([ini.zeta]);        
%     end

    iz = 1;
    M = OHC(iz).M;
    
    % Initial value of OHC membrane potential is given assuming
    % 0.45 resting open probability
    OHC(iz).Vnd0 = Vnd0((iz-1)*nv+1:iz*nv);
    OHC(iz).Vnd = OHC(iz).Vnd0;
    M.Vm0 = Vnd0((iz-1)*nv+2) - Vnd0((iz-1)*nv+3);
    M.Vm = M.Vm0;
    M = compute_membrane_capacitance(M, 0);

    OHC(iz).M = M;
    OHC(iz).ini.Gm = M.Gm;
    OHC(iz).ini.Cm = M.Cm;
    OHC(iz).ini.alpha = M.alpha;
    OHC(iz).ini.zeta = M.zeta;
    %% end of model OHC circuit

    %% initiate xHB

    ww = 2*pi*MP.freq;
    xx_t = 1e-3*amplitude*sin(ww*tt);
    xx_f = fourier_transform(xx_t);
    %%

    %% evaluate xa
    xa_f = (HB.kG*HB.gamma*xx_f)/(1i*ww/HB.kA + HB.kES);
    xa_t = inverse_fourier_transform(xa_f);
    %%

    %% evaluate po
    po_t = zeros(1,nt);

    for ti = 1:nt
        po_t(ti) = channel_2state(HB, xx_t(ti), xa_t(ti));
    end

    po_f = fourier_transform(po_t);

    %%

    %% evaluate voltage
    [Ce,Ge,Ie,V0] = mElectrical2_cmod_MSA(MP,OHC,IHC);
    
    voltage_f = zeros(nv,nt);
    voltage_f(:,ffsym == 0) = V0;
    voltage_t = inverse_fourier_transform(voltage_f);
    po0 = po_f(ffsym == 0); 
    while true
        voltage_t_n = voltage_t;
        GV_nonlin_t = zeros(nv,nt);
        for ti = 1:nt
            GV_nonlin_t(:,ti) = evaluate_GV_nonlin(MP,OHC,voltage_t(:,ti),po_t(ti),po0);
        end
        GV_nonlin_f = fourier_transform(GV_nonlin_t);
        idx = find(ffsym == MP.freq); % main frequency component
        voltage_f(:,idx) = (1i*ww*Ce + Ge)\(Ie + GV_nonlin_f(:,idx));
        voltage_f(:,nt-idx+2) = conj(voltage_f(:,idx));
        voltage_t = inverse_fourier_transform(voltage_f);
        err = norm(voltage_t(:) - voltage_t_n(:))/norm(voltage_t(:));
        if err < 1e-3
            break
        end
    end
end