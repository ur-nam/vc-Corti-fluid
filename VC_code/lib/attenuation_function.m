function w = attenuation_function(MP,zz,ndof,tdof,BC,b,m,za)
    
    if ~exist('b','var') || ~exist('m','var')
        b = MP.b;
        m = MP.m;
    end
%     zz = Nd.Z(:);
    f = zeros(length(zz),1);
    if ~exist('za','var')
        za = 11;
    end
    zz = zz*1e-3;
    zz = zz + abs(min(zz)); % bringing to positive longitude
    idx = zz > za;
    f(idx) = b*(zz(idx) - za).^m;
    w = 1./exp(-f);
    w = repmat(w,[1 ndof]);
    w = reshape(transpose(w),[],1);
    w = spdiags(w,0,tdof,tdof);
    dof = reshape(transpose(logical(BC)),[],1);
    w = w(dof,dof);

end