function U = inverse_fourier_transform(Uf)
% input Uf should correspond with wavenumbers, i.e. symmetric
    L = size(Uf,2);
    fft_coef = 2*ones(1,L); fft_coef(:,1) = 1;
    %shift back       
    Uf = [Uf(:,end/2+1:end),Uf(:,1:end/2)];
    U = real(ifft(Uf*L./fft_coef,[],2));
    %shift fft to correspond to wave numbers, more info in the personal notes   
end