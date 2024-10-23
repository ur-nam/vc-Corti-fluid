function MP = set_MP(z)

    global coef
    if isempty(z)

%         prompt{1} = 'Input directory';
%         prompt{2} = 'Output directory';
%         prompt{3} = 'Center location in [mm]';
%         prompt{4} = 'L of modeled piece in [mm]'; % DO NOT change this because this code uses pre-meshed fluid domain which assumes L = 4 mm.
%         prompt{5} = 'Stim freq range [f1, f2] in [kHz]';
%         prompt{6} = 'Num of frequencies [0: only input]';
%         prompt{7} = 'Stim pressure level in [dB SPL]';
%         prompt{8} = 'Permeation coef [Alpha, Beta]';        
%         prompt{9} = 'Stimulation protocol [0:Pure tone, 1:Multi tone]';
%         prompt{10} = 'Corti fluid [1:On, 0:Off]';
%         prompt{11} = 'Sensitive cochlea';
%         prompt{12} = 'Nonlinearity';
%         prompt{13} = 'Stim current level in [mA]';
%         
%         name='Input for Peaks function';
%         numlines=1;
        defaultanswer={'./hinput/','./houtput/','6','12',...
                '1','0','80','[1e3,1e0]','0','1','1','0','0.1'};
        app = InputPrompt(defaultanswer);
        waitfor(app,'output');
        z = app.output;
        delete(app);
% 
%         z = inputdlg(prompt,name,numlines,defaultanswer);
%         if isempty(z)
%             error('The run was terminated.');
%         end
    end
    
    for ii=3:length(z)
        z{ii} = str2num(z{ii}); % No str2double() won't work
    end


    idir = add_directory_divider(z{1});
    odir = add_directory_divider(z{2});
    MP.idir = idir;
    MP.odir = odir;
    
    if ~exist(odir,'dir') && ~isempty(odir)
        mkdir(odir);
    end

    MP.Estim = 0; % default: 0
    % basic parameters taken from the input %
    MP.loc = z{3};             % [mm]    
    if z{4} == 0
        MP.vMC = 1; %
        MP.loc = z{3} - 1; % 1mm offset to get PoI under slit % [mm]    
        MP.length_BM = 4000; % [um]
    else
        MP.vMC = 0; %
        MP.loc = z{3};             % [mm]   
        MP.length_BM = z{4}*1e3;   % [um]
    end
    MP.extension = 0;       % [um]
    MP.dZ = 10;
    MP.xx = MP.loc-0.5*MP.length_BM*1e-3:0.01:MP.loc+0.5*MP.length_BM*1e-3;
    MP.zz = (-0.5*(MP.length_BM + MP.extension):MP.dZ:0.5*(MP.length_BM + MP.extension))';
    MP.nz = length(MP.zz);
    MP.Nf = z{6};
    if MP.Nf == 0
        MP.freq = z{5};
        MP.Nf = length(MP.freq);
    else
        MP.freq = transpose(logspace(log10(z{5}(1)),log10(z{5}(2)),MP.Nf));
    end

    MP.BMC = {'A5'};
    MP.Pstim = z{7};            % dB SPL at the stapes re. 20 uPa
    if MP.vMC
        MP.Istim = z{13};
    end
    MP.opt_cfld = z{10}; % 1: consider corti fluid
    MP.p_alpha = coef.perm_multiplier*z{8}(1); % permeability constant k/(mu*L) [ms*um^2/ng]
    MP.p_beta = z{8}(2);
    MP.sensitive_cochlea = z{11};
    MP.nonlin = z{12}; % 1:nonlinear, 0:linear
    % hard-coded parameters %
    MP.H = 600;                % [um], This may be not used. Instead botH and topH are used.
    MP.CF_H = 50; % Cortfi fluid height [um]
    MP.rho = coef.rho_multiplier*1e-3; % fluid density [mg/mm^3]
    MP.nu = coef.nu_multiplier*700; % kinematic viscosity [um^2/ms]
    MP.m = 3; % power of attenuation function to resolve apical wave reflection
    MP.b = 1;

    % model options %
    MP.NQ = 4; % number of gauss points for fluid quadrilateral element
    MP.Visc = 0;        % Indicates viscous (1) or inviscid (0) model
    MP.key_FSI_2Chambers = 3;
    MP.key_FSI = 1;
    MP.opt_wtr_shnt = 0; % 1: consider water shunt
    MP.shnt_r_fact = 1; % water shunt relaxation factor
    MP.ind_MET = 2; % num of MET channel states, either 2 or 10
    MP.opt_in = 1;              % 1: pressure,  2: acceleration
    if z{9} == 1
        MP.opt_stim = 'multi'; % pure, multi
    elseif z{9} == 0
        MP.opt_stim = 'pure';
    end
    
    MP.phase = 0.5*pi;              % phase between M-stim and E-stim in degrees     
    MP.dtau = 1e0;%2.4142; % pseudo-time step

end

function pdir = add_directory_divider(pdir)
% % function pdir = add_directory_divider(pdir)
% % If there is no directory divider at the end of the string pdir
% % add it.

    if ~isempty(pdir),
        for ii=1:length(pdir),
            if strcmp(pdir(ii),'\'),
                pdir(ii) = '/';
            end        
        end

        if ~(strcmp(pdir(length(pdir)),'/') || strcmp(pdir(length(pdir)),'\')),
            pdir = [pdir '/'];
        end
    end
    
end