function set_Global_coefficients(z)
% function set_Global_coefficients(z)
% z: structure containing parameters; usually used for reruns
%
    
    clear global;
    global coef

    if exist('z','var')
        coef = z;
    else
        % gating swing
        coef.HB_b_A = 1.3;
        coef.HB_b_B = 1.3;
    
        %model_geometry
        coef.k_OHB_A = 1;
        coef.k_OHB_B = 1;
    
        coef.thick_A = 1;
        coef.thick_B = 1;
       
        coef.width_BM_A = 1.1;
        coef.width_BM_B = 1.0;
    
        coef.thick_PC_A = 1.0;
        coef.thick_PC_B = 1.0;
    
        coef.thick_BM_A = 0.95;
        coef.thick_BM_B = 1;
       
        coef.diam_DC_B = 1.0;
        coef.diam_DC_A = 1.4;

        coef.diam_OHC_B = 1.0;
        coef.diam_OHC_A = 1.2;

        coef.thick_TM_A = 1;
        coef.thick_TM_B = 1;
    
        %model_moduli
        coef.e_OHC = 1;
        coef.e_DCp = 1; % 0.1-1
        coef.e_DCb = 6.5;
        coef.e_OPC = 0.8;
        coef.e_IPC = 2;
        coef.e_TMx1_A = 1;
        coef.e_TMx1_B = 2;
        coef.e_RLx1 = 1;
        coef.e_RLp = coef.e_RLx1;
        coef.e_DCz = 1;
        coef.e_BMz = 1;
        coef.e_TMz = 1;
        coef.e_YMz = 1;
    
        %fluid
        coef.dcf_32_A = 0.6;
        coef.dcf_23_A = 0.6;
    
        coef.dcf_32_B = 0.6;
        coef.dcf_23_B = 0.6;
    
        %damping
        d_opt = 0;
        coef.beta_multiplier = 0.9;
        coef.nu_multiplier = 100; % Corti fluid nu
        if d_opt == 1
            coef.BM_damping_coef = 0.8;
            coef.RL_damping_coef = 1.0;
            coef.TM_damping_coef = 1.0;
            coef.DC_damping_coef = 1.2;
            coef.PC_damping_coef = 0.8;
        end
    
        %properties
        coef.rho_multiplier = 1; % fluid density
        coef.perm_multiplier = 1;

    end

        coef.hinge_ratio = 0.9;
        coef.OHC_damping = 1;
    
end