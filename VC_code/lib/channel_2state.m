function po = channel_2state(HB, xx, xa)

    xx = reshape(xx,1,[]);
    xa = reshape(xa,1,[]);

    dE = [HB.z].*(xx - xa - [HB.Xo]);
    kBT = 4e-3;
    kCO = [HB.kF].*exp(0.5*dE/kBT);      % rate from closed to open state
    kOC = [HB.kR].*exp(-0.5*dE/kBT);     % rate from open to closed state
    po = kCO./(kCO+kOC);

%     dT = fMP.dt;
%     if po + dpdt*dT >= 1
%         dpdt = 0.99*(1-p)/dT;
%     elseif po + dpdt*dT <= 0
%         dpdt = -0.99*p/dT;
%     end

end