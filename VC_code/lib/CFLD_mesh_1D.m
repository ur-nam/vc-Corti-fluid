function CFLD = CFLD_mesh_1D(CFLD,sNd)

    xx = sNd.Z(CFLD.nBot);
    Nd.width = CFLD.Nd.width;
    Nd.Z = xx;

    Nd.N = length(xx);
    El.node = [(1:Nd.N-1).',(2:Nd.N).'];
    El.N = size(El.node,1);

    idx = [find(strcmp(sNd.name,'BB')),find(strcmp(sNd.name,'DD')),...
           find(strcmp(sNd.name,'AM')),find(strcmp(sNd.name,'A5')),...
           find(strcmp(sNd.name,'AP')),find(strcmp(sNd.name,'A3')),...
           find(strcmp(sNd.name,'A2')),find(strcmp(sNd.name,'A1')),...
           find(strcmp(sNd.name,'A0')),find(strcmp(sNd.name,'C1')),...
           find(strcmp(sNd.name,'CC'))]; % entire OoC

%     idx = [find(strcmp(sNd.name,'BB')),find(strcmp(sNd.name,'DD')),...
%            find(strcmp(sNd.name,'AP')),find(strcmp(sNd.name,'A3')),...
%            find(strcmp(sNd.name,'A2')),find(strcmp(sNd.name,'A1')),...
%            find(strcmp(sNd.name,'A0')),find(strcmp(sNd.name,'C1')),...
%            find(strcmp(sNd.name,'CC'))]; % without DC

    idx = idx(CFLD.idx,:);

    Nd.p_num = size(idx,2);
    Nd.perim = idx;

    A = zeros(Nd.N,1);
    H = zeros(Nd.N,2*Nd.p_num);

    for ni = 1:Nd.N
        idx = Nd.perim(ni,:).';
        xx = [sNd.X(idx),sNd.Y(idx)];
            % formulation for approximate calculation of area change based on the shoelace formula
        A(ni,1) = polyarea(sNd.X(idx),sNd.Y(idx));
        H(ni,:) = -1*polyarea_change(xx); % the -1 is for the counterclock-wise direction of the node numbering, it could be reversed for convenience    
    end

    Nd.H = H;
    Nd.A = A;
    CFLD.Nd = Nd;
    CFLD.El = El;

end