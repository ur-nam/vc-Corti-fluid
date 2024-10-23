function [f,yy] = freq_sweep_plot

    epsilon = 1e-6;
    rdir = './results/freq_sweep_3/';
    list = ls(rdir);
    n = size(list,1);

    load('./input/OHC_ini.mat','OHC');
    icount = 0;
    for i = 1:n
        if ~contains(list(i,:),'NLFD'), continue; end
        load([rdir,list(i,:)],'MP');
        if ~((MP.stim - 30) < epsilon), continue; end        
        load([rdir,list(i,:)],'MP','Uf');
        icount = icount + 1;
        odof = MP.dof(1).n + (1:MP.dof(2).n);
        uf = Uf(odof,:);
        idx = (MP.ffsym == 0);
        uu = Uf(odof,idx);
        uu = reshape(uu,5,MP.nz);
        Vm0 = reshape(uu(2,:) - uu(3,:),MP.nz,1);
        
        % dVm time domain reconstruction
%         nt = 201;
%         tt = linspace(0,MP.tt(end),nt);
%         uu = real(uf(:,MP.ind)*exp(1i*2*pi*(MP.ffsym(MP.ind).').*tt));
%         uu = reshape(uu,5,MP.nz,nt);
%         Vm = reshape(uu(2,:,:) - uu(3,:,:),MP.nz,nt);
%         dVm = Vm - Vm0;
%         Vmax = max(dVm,[],2);

        % membrane voltage transfer function
        edof = MP.dof(1).edof;                
        idx = (MP.ffsym == MP.freq);
        po = Uf(edof,idx);
        num_section = [201,601,1001];
        for si = 1:length(num_section)
            n = num_section(si);
            voltage = reshape(uf(:,idx),5,MP.nz);
            dV = voltage(1,n) - voltage(2,n);
            I_MET(si,icount) = 3*(1i*(2*pi*MP.freq)*OHC(n).S.C + OHC(n).S.Gmax*po(n))*dV;
            dVm(si,icount) = voltage(3,n) - voltage(2,n);
            f_c(n) = 1/(2*pi*(OHC(n).M.Cm + OHC(n).M.alpha));
        end

        f(1,icount) = MP.freq;
        
    end

    yy = abs(dVm./I_MET);
    pp = angle(dVm./I_MET);

    figure;
    subplot(2,1,1)
    loglog(f,yy); hold on; xline(f_c);
    LPF = @(a,b,x) abs(a./(1 + x./1i*2*pi*b));
    c = fit(f.',yy(1,:).',LPF,'Start', [0.1,0.1]);
    loglog(f,feval(c,f),'k--');

    subplot(2,1,2)
    semilogx(f,unwrap(pp,[],2)/(2*pi)); hold on; xline(f_c);

%     figure;
%     subplot(2,1,1)
%     loglog(f,abs(I_MET)); hold on; xline(f_c);
%     subplot(2,1,2);
%     semilogx(f,unwrap(angle(I_MET),[],2)/(2*pi)); hold on; xline(f_c);
% 
%     figure;
%     subplot(2,1,1)
%     loglog(f,abs(dVm)); hold on; xline(f_c);
%     subplot(2,1,2);
%     semilogx(f,unwrap(angle(dVm),[],2)/(2*pi)); hold on; xline(f_c);
end