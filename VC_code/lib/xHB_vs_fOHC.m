function [xx,po,vm] = xHB_vs_fOHC(loc,freq)

    nx = 21;
    xx = linspace(0.1,10,nx);

    po = zeros(1,nx);
    vm = zeros(1,nx);
    for i = 1:nx
        [xx_t,xa_t,po_t,voltage_t,tt] = xHB_to_voltage(loc,freq,xx(i));
        po(i) = max(po_t);
        dVm = voltage_t(2,:)-voltage_t(3,:);
        vm(i) = max(dVm-mean(dVm));
    end

    figure; plot(xx,po);
    figure; plot(xx,vm);


end