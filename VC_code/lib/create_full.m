function Matrix = create_full(Nd,El,matrix)

    NDOF = 6;
    parts = zeros(144*El.N,3);
    cnt = 0; % counter variable

    for i = 1:El.N
    
        nd1 = El.Nd1(i);
        nd2 = El.Nd2(i);
        m = El.mat(i).(matrix);
        % for the beam elements
            if El.type(i) < 1

                rdof = [(nd1-1)*NDOF + (1:6),(nd2-1)*NDOF + (1:6)];
                
                % This function fills in the 3-column input for the sparse
                % function.
                [parts,cnt] = fillR(rdof,rdof,m,parts,cnt);
                
            else
                rdof = [(nd1-1)*NDOF + (1:3),(nd2-1)*NDOF + (1:3)];                
                [parts,cnt] = fillR(rdof,rdof,m,parts,cnt);
                
            end
    end

    % Remove excess zero elements resulting from 144 DoF assumption.
    zr = parts(:,1)==0;
    parts(zr,:) = [];
    % Use sparse function to create the matrix
    Matrix = sparse(parts(:,1),parts(:,2),parts(:,3),Nd.tdof,Nd.tdof);
end
