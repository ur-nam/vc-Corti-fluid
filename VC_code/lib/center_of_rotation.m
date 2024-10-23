function [radians,slenderness] = center_of_rotation(Nd,MP,R,P,C,loc)
% rotation offset of point P from center C

    NDOF = 6;
    P_node = strcmp(Nd.name,P);
    C_node = strcmp(Nd.name,C);
    udof = reshape(MP.dof(2).udof,NDOF,Nd.N);
    ux = R.Uf(udof(1,:),:);
    uy = R.Uf(udof(2,:),:);
    dx = ux(P_node,:);
    dy = uy(P_node,:);
    px = Nd.X(P_node);
    py = Nd.Y(P_node);
    cx = Nd.X(C_node);
    cy = Nd.Y(C_node);

    tt = linspace(0,2*pi,21);

    radians = zeros(1,MP.Nf);
    slenderness = zeros(1,MP.Nf);
    for fi = 1:MP.Nf
        xx = dx(:,fi);
        yy = dy(:,fi);

        if exist('loc','var')
            idx = knnsearch(MP.xx(:),loc);
        else           
            y_max = max(abs(yy));
            [~,idx] = min(abs(abs(yy) - 0.5*y_max));
        end

        dd = [real(xx(idx)*exp(1i*tt));real(yy(idx)*exp(1i*tt))];
        [d_max,ti] = max(vecnorm(dd,2,1));
        [d_min,~] = min(vecnorm(dd,2,1));
        slenderness(fi) = (d_max-d_min)/d_max;
        P_vec = transpose(dd(:,ti));
        PC_vec = [px(idx),py(idx)] - [cx(idx),cy(idx)];
        radians(fi) = asin(abs(P_vec(1)*PC_vec(2)-P_vec(2)*PC_vec(1))/norm(PC_vec)/norm(P_vec));
    end
end