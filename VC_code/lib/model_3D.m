function [Nd, El, MP] = model_3D(MP, dim, opt_lgradient, opt_MP_set)
% #
% # function [Nd, El, MP] = model_3D(MP, dim, opt_lgradient)
% #
% # Units are in pN, micrometer and Pa (N/m^2 or pN/um^2)
% #
% #
% # loc = location of partition 'b', 'a'
% # dim = 3, 2 (3dimension, 2 dimension)
% # opt_lgradient = optioin of longitudinal gradient 0 = no longitudinal
% #  gradient, 1 = gradient exists along longitudinal axis
% #
% # First written: 2009.06.xx
% # Last modified: 2021.07.xx
% #
% # by: Jong-Hoon Nam at University of Rochester
% #


MP.NDOF = 6;

if exist('dim','var')
    MP.dim = dim;   % 2D model or 3D model
else
    MP.dim = 3;
end


if isnumeric(MP.loc)
    if MP.loc <= 6
        location = 'b';
    else
        location = 'a';
    end
else
    fprintf(1,'\nError in model_3D, ''loc'' should be either numeric, ''a'' or ''b'' \n');
    MP.loc = 2;
    location = 'b';
end


if opt_MP_set == 0  
    MP = model_geometry(MP, location, opt_lgradient);    
    [Nd, MP] = model_nodes_og(MP);
    El = model_elements_og(Nd, MP);
    El.hinge = false(El.N,2);
    Nd.Master_Slave = zeros(sum(El.hinge(:)),2);
    % [Nd,El] = define_hinge_nodes(Nd,El);
    El = model_moduli5d_og(El, Nd, MP);    
else % This is a new set of geometry and material props updated on Apr 2021
    MP = model_geometry2021_test(MP, location, opt_lgradient);
    [Nd, MP] = model_nodes_test(MP);
    El = model_elements_test(Nd, MP);
    [Nd,El] = define_hinge_nodes(Nd,El);
    El = model_moduli5_2021_test(El, Nd, MP);    
end

% Nd = model_pArea(Nd, El);

MP.xOHC = 0;
MP.fOHC = 0;

end

function Nd = model_pArea(Nd, El)
% Compute the fluid-interacting surface area of nodes that are incorporated
% with membrane structures.
Nd.pArea = zeros(size(Nd.X));
plate_BM = sum(El.type>2.5)>0;

eE2 = strcmpi(El.name,'BMp')| strcmpi(El.name,'BMa') | strcmpi(El.name,'BMza') | strcmpi(El.name,'BMzp') ... % either BM
    | strcmpi(El.name,'TMx1') | strcmpi(El.name,'TMx2') ... % or TM
    | strcmpi(El.name,'TMz0') | strcmpi(El.name,'TMz1') ...
    | strcmpi(El.name,'RLx1') | strcmpi(El.name,'RLz');    % or RL


eE2p = strcmpi(El.name,'BMza') | strcmpi(El.name,'BMzp');
eE3p = strcmpi(El.name,'TMx1') | strcmpi(El.name,'TMx2') ... % TM
    | strcmpi(El.name,'TMz0') | strcmpi(El.name,'TMz1') ...
    | strcmpi(El.name,'RLx1') | strcmpi(El.name,'RLz'); % or RL



for ind = 1:Nd.N
    
    if plate_BM  % if the BMz is a plate element
        eE1 = El.Nd1==ind | El.Nd2==ind | El.Nd3==ind | El.Nd4==ind; % Element engaged to node ind
        
        
        if sum(eE2p)~=0
            eE = eE1 & eE2p;
            pA = (0.5*El.dim(eE,1)).*(0.5*El.dim(eE,2));    % pressurized area equals the length times the width
            Nd.pArea(ind) = sum(pA);
        elseif sum(eE3p)~=0
            eE = eE1 & eE3p;
            if sum(eE)==2, a=1/2; elseif sum(eE)==3, a=2/3;  else, a=1; end  %treat the boundaries
            pA = (0.5*El.L(eE)).*(0.5*El.dim(eE,1));    % pressurized area equals the length times the width
            Nd.pArea(ind) = a*sum(pA);
        else
        end
    else                  % if the BMz is a beam element
        eE1 = El.Nd1==ind | El.Nd2==ind; % Element engaged to the present node of interest (ind)
        eE = eE1 & eE2;
        if sum(eE)==2, a=1/2; elseif sum(eE)==3, a=2/3; else, a=1; end  %treat the boundaries
        pA = (0.5*El.L(eE)).*(0.5*El.dim(eE,1));    % pressurized area equals the length times the width
        Nd.pArea(ind) = a*sum(pA);
    end
    
end

end