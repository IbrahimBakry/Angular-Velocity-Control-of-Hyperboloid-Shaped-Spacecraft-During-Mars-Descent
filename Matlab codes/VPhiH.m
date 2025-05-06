function [out]=VPhiH(in)
    g = 9.81;
    q = ;
    s = ;
    m = ;
    cxv = ;
    
    v = in(1);
    theta = in(2);
    H = in(3);

    dvdt = -Cxv*q*s/m - g*sin(theta);
    dTdt = (-g*cos(theta))/v;
    dHdt = v*sin(theta);
    
    Cyn = Cyv*cos(alpha)-Cxv*sin(alpha);
    Cx = Cxv*cos(alpha)+Cyv*sin(alpha);
    Cyvalpha = dCyvdalpha;
    Mzn = mzn*q*s*L + q*s*L*(cos(phi)*(mzf+Cx*DeltaY) - sin(phi)*(myf-Cx*DeltaZ))
    Myn = q*s*L*(cos(phi)*(myf-Cx*DeltaZ) - sin(phi)*(mzf+Cx*DeltaY))
    Mx = q*s*L*Cyn*(DeltaZ*cos(phi) + DeltaY*sin(phi))
    Mxv = Mx*cos(alpha) - Myn*sin(alpha);
    
    Mznphi = -Cx1*DeltaZn*q*s + (mzf*cos(phin)-myf*sin(phin));
    N = mu*Mzn*dqdt/(I*q)+mu*Mznphi*dphidt/I + mu*(Ix*OmegaX-2*OmegaG*cos(alpha))*Myn/I - mu*Mx*OmegaG*sin(alpha)/I - mu*Cyv*q*s*(Ix*OmegaX-OmegaG*cos(alpha))*(Ix*OmegaX-2*OmegaG*cos(alpha))/(m*v);
    D = (OmegaG^2 + (Ix*OmegaX-OmegaG*cos(alpha))*(Ix*OmegaX-2*OmegaG*cos(alpha)));
    dadt = N./D;
    dOmegaxdt = mu*Mx/Ix;
    dPhidt = OmegaX - OmegaY*cos(alpha);
    dGammadt = (Ix*OmegaX/(2*cos(alpha))) + (Ix^2*OmegaX^2/4 - Mzn*cotd(alpha)/I)^0.5/cos(alpha);
    
    out = [dvdt; dTdt; dHdt];
end