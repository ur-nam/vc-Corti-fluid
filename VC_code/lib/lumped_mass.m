function  El = lumped_mass(Nd, El,  JH)
% The masses of BM and OC are lumped to the elements along the z-axis, BMz and DCz

    
% %     mthick_BMz = thickness_BMz(Nd, El, eBMz, JH);
    
    eBMza = strcmp('BMza',El.name);
    eBMzp = strcmp('BMzp',El.name);
    eBMz = eBMza | eBMzp;
    mthick_BMza = thickness_BMz(Nd, El, eBMza, JH);
    mthick_BMzp = thickness_BMz(Nd, El, eBMzp, JH);
    
    eRLz = strcmp('RLz',El.name);
    eX = eRLz;
    mthick_RLz = thickness_Z(Nd, El, eX, 'thick_BMz', JH);
    mfactor = 0.2*mthick_RLz./El.dim(eX,2); % This is not good. This should be defined in model_geometry and in JH.    
    El.rho(eRLz) = El.rho(eRLz) + 0*mfactor*1e-3;
    
    eTMz0 = strcmpi(El.name,'TMz0');
    mthick_TMz0 = thickness_Z(Nd, El, eTMz0, 'thick_TMza', JH);
    
    eTMz1 = strcmpi(El.name,'TMz1');
    mthick_TMz1 = thickness_Z(Nd, El, eTMz1, 'thick_TMzb', JH);    
   
    if mean(El.type(eBMz))>2.5,
        
        eX = eBMz; mthick = mthick_BMz;
        mfactor = mthick./El.thick(eBMz);        
        El.rho(eX) = mfactor*1e-3;
    
    else        
        eX = eBMzp; mthick = mthick_BMzp;
        
        
        
        
        mfactor = 1.5*mthick./El.dim(eX,2); % factor 1.5 in order to include supporting cells in OC
        
        
        
        
        El.rho(eX) = El.rho(eX) + mfactor*1e-3; % this needs a revisit
        
        eX = eBMza; mthick = mthick_BMza;
        mfactor = mthick./El.dim(eX,2);
       
        El.rho(eX) = El.rho(eX) + 0*mfactor*1e-3;
        
    end    
    
    eX = eTMz0; mthick = mthick_TMz0;
    mfactor = mthick./El.dim(eX,2);
    El.rho(eX) = mfactor*1e-3;

    eX = eTMz1; mthick = mthick_TMz1;
    mfactor = mthick./El.dim(eX,2);
    El.rho(eX) = mfactor*1e-3;
    
end


function thick_BMz = thickness_BMz(Nd, El, eX, JH)
    
    qJH = consider_geometrical_lgradient(JH);
    
    if mean(El.type(eX))>2.5,   % if it is a plate
        nds = [El.Nd1(eX), El.Nd2(eX), El.Nd3(eX), El.Nd4(eX)];
    else                        % else if it is a beam
        nds = [El.Nd1(eX), El.Nd2(eX)];
    end
    xi = mean(Nd.X(nds),2);
    zi = mean(Nd.Z(nds),2);
    if isfield(JH,'zz')
        zz = JH.zz;
    else
        zz = (-0.5*JH.length_BM:JH.dZ:0.5*JH.length_BM)';
    end
    peak_thick = interp1(zz,qJH.thick_BMz,zi);
%     thick_fluid = interp1(zz,qJH.thick_fluid,zi);

    thick_fluid = interp1(zz,qJH.thick_fluid,0*zi);
    % fluid-mass gradient may be ignored (at the singe stimulating frequency, the fluid mass may not increase along the length)
    
    
    e_BMwidth = interp1(zz,qJH.width_BM,zi);
        
    xi = xi./e_BMwidth;    
    
    x_profile = zeros(size(xi));
    za = xi<1/3; % zona arcuate
    zp = xi>1/3; % zona pectinate
    ze = xi==0 | xi==1; % edge
    x_profile(za) = -36*xi(za).*(xi(za)-1/3); 
    x_profile(zp) = -9*(xi(zp)-1/3).*(xi(zp)-1);
    thick_BMz = x_profile.*peak_thick;
    thick_BMz(ze) = 0.1;
    
    thick_BMz = thick_BMz + 0.5*thick_fluid;
    
end



function thick_TMz = thickness_Z(Nd, El, eX, prop, JH)
    
    qJH = consider_geometrical_lgradient(JH);
    
    if mean(El.type(eX))>2.5,   % if it is a plate
        nds = [El.Nd1(eX), El.Nd2(eX), El.Nd3(eX), El.Nd4(eX)];
    else                        % else if it is a beam
        nds = [El.Nd1(eX), El.Nd2(eX)];
    end
    
    zi = mean(Nd.Z(nds),2);

    if isfield(JH,'zz')
        zz = JH.zz;
    else
        zz = (-0.5*JH.length_BM:JH.dZ:0.5*JH.length_BM)';
    end
    
    zz = zz + mean([Nd.Z(El.Nd1(eX)); Nd.Z(El.Nd2(eX))]);
    
% %     thick_fluid = interp1(zz,qJH.thick_fluid,zi);    
    
    thick_fluid = interp1(zz,qJH.thick_fluid,0*zi);
    % fluid-mass gradient may be ignored (at the single stimulating frequency, the fluid mass may not increase along the length)
    
    thick_TMz = interp1(zz,qJH.(prop),zi);    
    
    thick_TMz = thick_TMz + 0.5*thick_fluid;

end
