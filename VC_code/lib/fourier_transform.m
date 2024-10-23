function Uf = fourier_transform(U)
    N = ndims(U);
    L = size(U,N);
    fft_coef = 2*ones(size(U)); 
    if N == 2
        fft_coef(:,1) = 1;
    elseif N == 3
        fft_coef(:,:,1) = 1;
    end
    Uf = fft_coef.*fft(U,[],N)/L;
    %shift fft to correspond to wave numbers
    if N == 2
        Uf = [Uf(:,end/2+1:end),Uf(:,1:end/2)];
    elseif N ==3
        Uf = cat(3,Uf(:,:,end/2+1:end),Uf(:,:,1:end/2));
    end        
end