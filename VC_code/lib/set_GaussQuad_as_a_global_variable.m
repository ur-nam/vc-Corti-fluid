function [xi, eta, w] = set_GaussQuad_as_a_global_variable(m)

    global GaussQuad;
    
    %-------
    if(m==4)
    %-------
    
    al = 0.577350269189626;
    xi(1) = -al; eta(1) = -al; w(1) = 1;
    xi(2) = -al; eta(2) = al; w(2) =  1;
    xi(3) = al; eta(3) = -al; w(3) =  1;
    xi(4) = al; eta(4) = al; w(4) =  1;

    %-----------
    elseif(m==9)
    %-----------

    al = 0.774596669241483;
    o1 = 0.555555555555556;
    o2 = 0.888888888888889;

    xi(1) = -al;
    xi(2) = -al;
    xi(3) = -al;
    xi(4) = 0;
    xi(5) = 0;
    xi(6) = 0;
    xi(7) = al;
    xi(8) = al;
    xi(9) = al;

    eta(1) = -al;
    eta(2) = 0;
    eta(3) = al;
    eta(4) = -al;
    eta(5) = 0;
    eta(6) = al;
    eta(7) = -al;
    eta(8) = 0;
    eta(9) = al;

    w(1) = o1*o1;
    w(2) = o1*o2;
    w(3) = o1*o1;
    w(4) = o1*o2;
    w(5) = o2*o2;
    w(6) = o1*o2;
    w(7) = o1*o1;
    w(8) = o1*o2;
    w(9) = o1*o1;

    %--
    end
    %--

    %-----
    % done
    %-----
    GaussQuad.xi = xi;
    GaussQuad.eta = eta;
    GaussQuad.w = w;
end
