function [nF, fdof] = make_unit_force2(Nd, fnode, fdir)
%      [nF, fdof] = make_unit_force(Nd, fnode, fdir)
% 
%      Program unit: generate the force vector
%      Description: 
%      Programmed: 2009.10.13
%      Last Modified: 
%      by Jong-Hoon Nam
% -----------------------------------------------------------------------
    mapF2R = Nd.mapF2R;

%     mapF2R = find_mapF2R(Nd);

    NDOF = 6;
    nF=zeros(max(mapF2R),1);
    
    if length(fdir) == 1,   % fdir indicates direction, 1=dx, 4=rx
        fdof = mapF2R(NDOF*(fnode-1) + fdir);
        if fdof~=0,
            nF(fdof) = 1; 
        end
    else                    % fdir indicates direction vector [cx, cy, cz]

        fdof = mapF2R(NDOF*(fnode-1) + (1:3))';
        temp = find(fdof~=0);
        nF(fdof(temp)) = fdir(temp);

% %         fdof = mapF2R(NDOF*(fnode-1)) + (1:3)';
% %         nF(fdof) = fdir;
    end
end
    
        
  