function dx = HB_deflection(MP,El,Uf,threshold,icount)

    ndof = 6;
    HBdir = 1;

    nANK = find(strcmp(El.name,'ANK'));
    nd1 = El.Nd1(nANK); nd2 = El.Nd2(nANK);
%     nOHB = find(strcmpi('OHB',El.name));
%     nd1 = El.Nd1(nOHB); nd2 = El.Nd2(nOHB);

    udof = MP.dof(2).udof;

    uu_f = Uf(udof,:);

    dx1_f = uu_f((nd1-1)*ndof + HBdir,:);
    dx2_f = uu_f((nd2-1)*ndof + HBdir,:);

    dx_f = dx2_f - dx1_f;
    dx = inverse_fourier_transform(dx_f);

%     % amplify initial condition
%     if icount == 2
%         a = 0.6*(100 - max(MP.Pstim)); % dB
%         dx = dx*(10^(a/20));
%     end

    % filter the relatively small displacements
    wave_envelope = max(dx,[],2);
    ref = max(wave_envelope);
    idx = 20*log10(wave_envelope/ref) + threshold > 0;
    dx(~idx,:) = 0;
end