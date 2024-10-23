function El = model_elements_test(Nd, MP)
% # function El = model_elements(Nd, MP)
% #
% # Define elements
% #
% # First written: 2009.06.xx
% # Last modified: 2013.11.04
% #
% # Written by: Jong-Hoon Nam at University of Rochester
% #
    
    qMP = consider_geometrical_lgradient(MP);

    nz = size(Nd.Z,1);
        
    El.Nd1=[]; El.Nd2=[]; El.Nd3=[]; El.Nd4=[]; El.name=[];
    

    nBM = MP.NX;

    kAP = find(strcmpi(Nd.name(1,:),'AP')); % x-index of the node at OPC root
    kAM = find(strcmpi(Nd.name(1,:),'AM')); % x-index of the node at DC root

    ix = 1:nBM;
    nd1 = (1:nz)' + (ix-1)*nz;
    nd2 = (1:nz)' + (ix)*nz;

    El.Nd1 = [El.Nd1; nd1(:)];
    El.Nd2 = [El.Nd2; nd2(:)];
    xOPC = qMP.OPC_root_location(1:nz).*qMP.width_BM(1:nz)*ones(1,nBM);
    xE = mean(Nd.X([nd1(:), nd2(:)]),2);
    name = repmat({'BMp'},(nz)*(nBM),1);
    name(xE<=xOPC(:)) = {'BMa'};
    El.name = [El.name;name];

    for iz = 1:nz  % section number along the longitudinal directon (z-dir)
    
        % BM: basilar membrane
        % IPC: inner pillar cell
        % OPC: outer pillar cell
        % DC:  Deiter's cell
        % RL: reticular laminar
        % OHC: outer hair cell
        % OHB: outer hair cell bundle
        % TM: tectorial membrane
        
        % % Basilar Membrane in radial direction
        
%         nd1 = ((1:nBM)'-1)*nz+ iz;
%         nd2 = ((1:nBM)')*nz+ iz;
%         
%         El.Nd1 = [El.Nd1; nd1];
%         El.Nd2 = [El.Nd2; nd2];
%         
%         for jj = 1:length(nd1)
%             if mean(Nd.X([nd1(jj), nd2(jj)])) < qMP.width_BM(iz)*(3/10) 
%                 El.name = [El.name;{'BMa'}];        % Basilar membrane acuate zone
%             else
%                 El.name = [El.name;{'BMp'}];        % Basilar membrane pectinate zone
%             end
%         end

        a1 = find(strcmp(Nd.name(1,:),'A1'));
        a3 = find(strcmp(Nd.name(1,:),'A3'));
        b1 = find(strcmp(Nd.name(1,:),'B1'));
        b2 = find(strcmp(Nd.name(1,:),'B2'));
        c1 = find(strcmp(Nd.name(1,:),'C1'));   
        c2 = find(strcmp(Nd.name(1,:),'C2'));
        c3 = find(strcmp(Nd.name(1,:),'C3'));
        c4 = find(strcmp(Nd.name(1,:),'C4'));
        c5 = find(strcmp(Nd.name(1,:),'C5'));
        c6 = find(strcmp(Nd.name(1,:),'C6'));
        dd = find(strcmp(Nd.name(1,:),'DD'));
        d1 = find(strcmp(Nd.name(1,:),'D1'));
        cc = find(strcmp(Nd.name(1,:),'CC'));        
        bb = find(strcmp(Nd.name(1,:),'BB'));
        ff = find(strcmp(Nd.name(1,:),'FF'));
        e1 = find(strcmp(Nd.name(1,:),'E1'));
        e2 = find(strcmp(Nd.name(1,:),'E2'));
        e3 = find(strcmp(Nd.name(1,:),'E3'));
        e4 = find(strcmp(Nd.name(1,:),'E4'));
        e5 = find(strcmp(Nd.name(1,:),'E5'));
        e6 = find(strcmp(Nd.name(1,:),'E0'));        
        pp = find(strcmp(Nd.name(1,:),'PP'));   
        
        a1 = (a1-1)*nz + iz; a3 = (a3-1)*nz + iz;
        b1 = (b1-1)*nz + iz; b2 = (b2-1)*nz + iz;
        c1 = (c1-1)*nz + iz; c2 = (c2-1)*nz + iz; c3 = (c3-1)*nz + iz; 
        c4 = (c4-1)*nz + iz; c5 = (c5-1)*nz + iz; c6 = (c6-1)*nz + iz; 
        dd = (dd-1)*nz + iz;
        d1 = (d1-1)*nz + iz;
        cc = (cc-1)*nz + iz; bb = (bb-1)*nz + iz; pp = (pp-1)*nz + iz;
        ff = (ff-1)*nz + iz;
        e1 = (e1-1)*nz + iz; e2 = (e2-1)*nz + iz; e3 = (e3-1)*nz + iz;
        e4 = (e4-1)*nz + iz; e5 = (e5-1)*nz + iz; e6 = (e6-1)*nz + iz;

        
        
%         bbp = bb + MP.tilt_OHC; % bbp: tilt of DCp toward the apical direction by distance of tilt_OHC
        
        bbp = pp + MP.tilt_OHC; % bbp: tilt of DCp toward the apical direction by distance of tilt_OHC
        bbq = pp - MP.tilt_OHC; % bbp: tilt of DCp toward the basal direction by distance of tilt_OHC
        
%         bbp = bb + MP.tilt_OHC; % bbp: tilt of DCp toward the apical direction by distance of tilt_OHC
%         bbq = bb - MP.tilt_OHC; % bbp: tilt of DCp toward the basal direction by distance of tilt_OHC
        
        a0 = 0 +  iz; % root of IPC
        
        ap = (kAP-1)*nz + iz; % root of OPC
        a2 = (kAM-1)*nz + iz; % root of DC
                        
        if rem(abs(Nd.Z(iz,1)),10)<0.1
            % % Pillar Cells
            if ~isempty(c1)
                El.Nd1 = [El.Nd1; a0; c1; c2; c3; ap; c4; c5; c6];
                El.Nd2 = [El.Nd2; c1; c2; c3; cc; c4; c5; c6; cc];
                El.name = [El.name; {'IPC'}; {'IPC'}; {'IPC'}; {'IPC'};...
                    {'OPC'}; {'OPC'}; {'OPC'}; {'OPC'};];
                %brackets
                El.Nd1 = [El.Nd1; c1; c3; a3; c6];
                El.Nd2 = [El.Nd2; a1; b1; c4; b1];
                El.name = [El.name; {'IPC'};{'IPC'};{'OPC'};{'OPC'};];
            else
                El.Nd1 = [El.Nd1; a0; ap];
                El.Nd2 = [El.Nd2; cc; cc];
                El.name = [El.name; {'IPC'}; {'OPC'};];                
            end

            if MP.tilt_OHC>=0 && nz-iz >= MP.tilt_OHC
                El.Nd1 = [El.Nd1; a2; dd; dd; bb];
                El.Nd2 = [El.Nd2; dd; bbp; bb; pp];
                El.name = [El.name; {'DCb'}; {'DCp'};{'OHC'};{'RLp'}];         
            elseif MP.tilt_OHC<=0 && iz >= -MP.tilt_OHC                
                El.Nd1 = [El.Nd1; a2; dd; dd; bb];
                El.Nd2 = [El.Nd2; dd; bbp; bb; pp];
                El.name = [El.name; {'DCb'}; {'DCp'};{'OHC'};{'RLp'}]; 
            else
                El.Nd1 = [El.Nd1; a2; dd;];
                El.Nd2 = [El.Nd2; dd; bb;];
                El.name = [El.name; {'DCb'}; {'OHC'}];
            end            
                        
            if MP.tilt_OHC>=0 && nz-iz < MP.tilt_OHC   % end of triangular repeats, suppressed
                El.Nd1 = [El.Nd1; dd];
                El.Nd2 = [El.Nd2; bbq];
                El.name = [El.name; {'DCp'};]; % end of triangular repeats, suppressed
            elseif  MP.tilt_OHC<=0 && iz < MP.tilt_OHC
                El.Nd1 = [El.Nd1; dd;];
                El.Nd2 = [El.Nd2; bbq;];
                El.name = [El.name; {'DCp'}];            
            end

            % % Reticular lamina            
            El.Nd1 = [El.Nd1; cc; b1];
            El.Nd2 = [El.Nd2; b1; bb];
            El.name = [El.name; {'RLx1'}; {'RLx1'}];


            % % OHC bundle and linear spring
            El.Nd1 = [El.Nd1; bb;];
            El.Nd2 = [El.Nd2; e1;];
            El.name = [El.name; {'OHB'};];

            % % Pseudo rod at OHC bundle location
            El.Nd1 = [El.Nd1; bb;];
            El.Nd2 = [El.Nd2; ff;];
            El.name = [El.name; {'xHB'};];
            
            % % Pseudo link to conveniantly calculate translational displ. 
            % % or damping of HB            
            El.Nd1 = [El.Nd1; ff;];
            El.Nd2 = [El.Nd2; e1;];
            El.name = [El.name; {'ANK'};];            
        end
        
        % % Tectorial Membrane
        El.Nd1 = [El.Nd1; e1; e2; e3; e4; e5;];
        El.Nd2 = [El.Nd2; e2; e3; e4; e5; e6;];
        
        
        El.name = [El.name; {'TMx2'}; {'TMx2'}; {'TMx2'}; {'TMx2'}; {'TMx1'};]; % 2/3 body zone, 1/3 attachement zoen 
% %         El.name = [El.name; {'TMx2'}; {'TMx2'}; {'TMx2'}; {'TMx1'}; {'TMx1'};]; % 1/2 body zone, 1/2 attachement zone
        % % radial stiffner along OHC basal ends
% %         El.Nd1 = [El.Nd1; c2];
% %         El.Nd2 = [El.Nd2; dd];
% %         El.name = [El.name;{'DCr'}];

        
    end
    
   
    if nz>1
        % % Tectorial Membrane in longitudinal direction
                     
        ee = find(strcmp(Nd.name(1,:),'E1'));
        nd1 = (ee-1)*nz + (1:nz-1)';     % Through E1
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1];
        El.Nd2 = [El.Nd2; nd2];

        El.name = [El.name;repmat({'TMz0'},nz-1,1)];
%         for iz=1:nz-1
%             El.name = [El.name;{'TMz0'}];
%         end

        ee = find(strcmp(Nd.name(1,:),'E2'));
        nd1 = (ee-1)*nz + (1:nz-1)';      % Through E2
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1];
        El.Nd2 = [El.Nd2; nd2];

        El.name = [El.name;repmat({'TMz1'},nz-1,1)];
%         for iz=1:nz-1
%             El.name = [El.name;{'TMz1'}];
%         end

        ee = find(strcmp(Nd.name(1,:),'E3'));
        nd1 = (ee-1)*nz + (1:nz-1)';      % Through E3
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1];
        El.Nd2 = [El.Nd2; nd2];

        El.name = [El.name;repmat({'TMz1'},nz-1,1)];
%         for iz=1:nz-1
%             El.name = [El.name;{'TMz1'}];
%         end    

        ee = find(strcmp(Nd.name(1,:),'E4'));
        nd1 = (ee-1)*nz + (1:nz-1)';      % Through E4
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1];
        El.Nd2 = [El.Nd2; nd2];

        El.name = [El.name;repmat({'TMz1'},nz-1,1)];
%         for iz=1:nz-1
%             El.name = [El.name;{'TMz1'}];
%         end

% %         ee = find(strcmp(Nd.name(1,:),'E5'));
% %         nd1 = (ee-1)*nz + (1:nz-1)';      % Through E5
% %         nd2 = nd1 + 1;
% %         El.Nd1 = [El.Nd1; nd1];
% %         El.Nd2 = [El.Nd2; nd2];
% % 
% %         for iz=1:nz-1,
% %             El.name = [El.name;{'TMz2'}];
% %         end    
        

        % % longitudinal stiffner along TC summits
        ee = find(strcmp(Nd.name(1,:),'CC'));
        nd1 = (ee-1)*nz + (1:nz-1)';
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1];
        El.Nd2 = [El.Nd2; nd2];

        El.name = [El.name;repmat({'TCz'},nz-1,1)];
%         for iz=1:nz-1
%             El.name = [El.name;{'TCz'}];
%         end 

        % % longitudinal stiffner along RL at hair bundle roots
        ee = find(strcmp(Nd.name(1,:),'BB'));
        nd1 = (ee-1)*nz + (1:nz-1)';
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1];
        El.Nd2 = [El.Nd2; nd2];

        El.name = [El.name;repmat({'RLz'},nz-1,1)];

        ee = find(strcmp(Nd.name(1,:),'B1'));
        nd1 = (ee-1)*nz + (1:nz-1)';
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1];
        El.Nd2 = [El.Nd2; nd2];

        El.name = [El.name;repmat({'RLz'},nz-1,1)];
%         for iz=1:nz-1
%             El.name = [El.name;{'RLz'}];
%         end
        
        ee = find(strcmp(Nd.name(1,:),'PP'));
        nd1 = (ee-1)*nz + (1:nz-1)';
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1];
        El.Nd2 = [El.Nd2; nd2];

        El.name = [El.name;repmat({'RLpz'},nz-1,1)];
%         for iz=1:nz-1
%             El.name = [El.name;{'RLpz'}];
%         end
        
        % % longitudinal stiffner along OHC basal ends
        ee = find(strcmp(Nd.name(1,:),'DD'));
        nd1 = (ee-1)*nz + (1:nz-1)';
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1];
        El.Nd2 = [El.Nd2; nd2];

        El.name = [El.name;repmat({'DCz'},nz-1,1)];
%         for iz=1:nz-1
%             El.name = [El.name;{'DCz'}];
%         end
        
    end % if nz>1,
    
    

    % % Basilar Membrane in longitudinal direction
    if ~isfield(MP,'opt_BM')
        fprintf(1,'\nBM element type unspecified. The program will use default beam element.\n');
        warndlg({'BM element type unspecified'; 'The program will use default beam element.'},'Never mind');
        MP.opt_BM = 0.1;
    end
       
    if MP.opt_BM>=3 % Four-node plate element will be used for the BM element in the z-direction
    
        El.Nd3 = zeros(size(El.Nd1));    
        El.Nd4 = zeros(size(El.Nd1)); 
      
        for ix = 1:nBM
            nd1 = (1:nz-1)' + (ix-1)*nz;
            nd2 = nd1 + 1;
            nd4 = (1:nz-1)' + (ix-0)*nz;
            nd3 = nd4 + 1;

            El.Nd1 = [El.Nd1; nd1];
            El.Nd2 = [El.Nd2; nd2];
            El.Nd3 = [El.Nd3; nd3];
            El.Nd4 = [El.Nd4; nd4];

            for iz=1:length(nd1)
                
                xOPC = qMP.OPC_root_location(iz)*qMP.width_BM(iz);                
                xE = mean(Nd.X([nd1(iz), nd2(iz), nd3(iz), nd4(iz)]));
                
                if xE<=xOPC
                    El.name = [El.name;{'BMza'}];
                else
                    El.name = [El.name;{'BMzp'}];
                end
                
            end
        end
    else    % Two-node beam element will be used for the BM element in the z-direction
        ix = 1:nBM+1;
        nd1 = (1:nz-1)' + (ix-1)*nz;
        nd2 = nd1 + 1;
        El.Nd1 = [El.Nd1; nd1(:)];
        El.Nd2 = [El.Nd2; nd2(:)];
        xOPC = qMP.OPC_root_location(1:nz-1).*qMP.width_BM(1:nz-1)*ones(1,nBM+1);
        xE = mean(Nd.X([nd1(:), nd2(:)]),2);
        name = repmat({'BMzp'},(nz-1)*(nBM+1),1);
        name(xE<=xOPC(:)) = {'BMza'};
        El.name = [El.name;name];
%         for ix = 1:nBM+1
%             nd1 = (1:nz-1)' + (ix-1)*nz;
%             nd2 = nd1 + 1;
% 
%             El.Nd1 = [El.Nd1; nd1];
%             El.Nd2 = [El.Nd2; nd2];
% 
%             for iz=1:length(nd1)
%                 
%                 xOPC = qMP.OPC_root_location(iz)*qMP.width_BM(iz);
%                 xE = mean(Nd.X([nd1(iz), nd2(iz)]));
%                 
%                 if xE<=xOPC
%                     El.name = [El.name;{'BMza'}];
%                 else
%                     El.name = [El.name;{'BMzp'}];
%                 end
%                 
%             end        
%         end                
    end
    
    El.N = length(El.Nd1);    
    
end
    