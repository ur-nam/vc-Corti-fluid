function convection = evaluate_CFLD_nonlinear_terms...
    (FLD, CU, CV)

    CFLD = FLD.CFLD;
    El = CFLD.El;
    convection = zeros(El.N,2*El.NQ);

    for kk = 1:El.N

        knodes = El.node(kk,:);
        convection(kk,:) = element_CFLD_conv(kk, El, CU(knodes), CV(knodes));
    end

end

function CCONV = element_CFLD_conv(ie, El, CU, CV)

    NQ = El.NQ;
    CCONV = zeros(1,2*NQ);

    for iq = 1:NQ
        cpsi = El.elm_psi(:,iq,ie);
        cgpsix = El.elm_gpsi(1:El.type,iq,ie);
        cgpsiy = El.elm_gpsi(1 + El.type:end,iq,ie);       
        
        u0 = cpsi'*CU(:);
        v0 = cpsi'*CV(:);        

        du0dx = cgpsix'*CU(:);
        du0dy = cgpsiy'*CU(:);
        dv0dx = cgpsix'*CV(:);
        dv0dy = cgpsiy'*CV(:);
        
        adv_x = u0*du0dx + v0*du0dy;
        adv_y = u0*dv0dx + v0*dv0dy;
        CCONV(1,iq) = adv_x; CCONV(1,iq + NQ) = adv_y;        
    end

end