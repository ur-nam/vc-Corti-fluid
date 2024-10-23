function res_val = evaluate_err(fMP,MP,Nd,PoI,Xf_n,Xf,icount,h)

    epsi = 1e-9;

    pdof = MP.dof(1).pdof;
%     udof = MP.dof(2).udof;
    edof = MP.dof(3).edof;
    adof = MP.dof(4).adof;
    odof = MP.dof(5).odof;    

    % time domain
    ndof = 6;
    udof = reshape(MP.dof(2).udof,ndof,Nd.N);
    udof_yBM = udof(2,PoI.zBM);
    yBM = inverse_fourier_transform(Xf(udof_yBM,:));
    yBM_n = inverse_fourier_transform(Xf_n(udof_yBM,:));
    c = max(yBM,[],2);
    r = max(yBM_n,[],2);   
    [~,idx] = max(c); 
    x_res_val = norm(c(:) - r(:))/(epsi + norm(c(:)));
    c = inverse_fourier_transform(Xf(odof,:));
    r = inverse_fourier_transform(Xf_n(odof,:));
    po_res_val = norm(c(idx,:) - r(idx,:))/(epsi + norm(c(idx,:)));

    % frequency domain
    n = size(fMP.ind,2); 
    c = abs(Xf(udof_yBM,fMP.ind(1:n)));
    r = abs(Xf_n(udof_yBM,fMP.ind));

%     % two-tone suppression
%     if (icount > 2) && exist('h','var')
%         fi = knnsearch(fMP.ffsym(fMP.ind).',MP.freq(2));
%         [~,idx] = max(c(:,fi));
%         yBM = abs(Xf(udof_yBM,:));
%         % probe f2
%         ff0 = knnsearch(fMP.ffsym.',MP.freq(2));
%         % f2-f1
%         ff1 = knnsearch(fMP.ffsym.',(MP.freq(2)-MP.freq(1)));
%         % f2+f1
%         ff2 = knnsearch(fMP.ffsym.',(MP.freq(2)+MP.freq(1)));
%         % 2f2-f1
%         ff3 = knnsearch(fMP.ffsym.',(2*MP.freq(2)-MP.freq(1)));
%         addpoints(h(3),icount,20*log10(yBM(idx,ff1)/yBM(idx,ff0)));
%         addpoints(h(4),icount,20*log10(yBM(idx,ff2)/yBM(idx,ff0)));
%         addpoints(h(5),icount,20*log10(yBM(idx,ff3)/yBM(idx,ff0)));
%         drawnow;        
%     end

    x_f_res_val = norm(c(:)-r(:))/(epsi + norm(c(:)));
    c = abs(Xf(odof,fMP.ind(1:n))); 
    r = abs(Xf_n(odof,fMP.ind));
    po_f_res_val = norm(c(:) - r(:))/(epsi + norm(c(:)));

    fprintf(1,'     yBM: %.3e, p_o: %.3e\n',x_f_res_val,po_f_res_val);
    if (icount > 2) && exist('h','var')% && x_f_res_val < 0.5)
        addpoints(h(1),icount,x_f_res_val);
        addpoints(h(2),icount,po_f_res_val);
        drawnow;
    end
    res_val = max(x_f_res_val,po_f_res_val);

end