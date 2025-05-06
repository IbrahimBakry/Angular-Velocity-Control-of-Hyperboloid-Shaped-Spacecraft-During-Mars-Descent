%% differentiation
syms mzn q s L I alpha Ixd omegax

omega = sqrt(-mzn*q*s*L*cot(alpha)/I);
omegaa = sqrt(Ixd^2*omegax^2/4 + omega^2);

omega12 = Ixd*omegax/2 + omegaa; 
% omega12 = Ixd*omegax/2 - omegaa;

D = omegax - omega12; 


%% dbdomega
syms mxa1 mxa2 ma1 ma2 Fa
dDdomegax = 1 - omega12;
dDdalpha = -(L*mzn*q*s*(cot(alpha)^2 + 1))/(2*I*((Ixd^2*omegax^2)/4 - (L*mzn*q*s*cot(alpha))/I)^(1/2));
% dDdalpha = (L*mzn*q*s*(cot(alpha)^2 + 1))/(2*I*((Ixd^2*omegax^2)/4 - (L*mzn*q*s*cot(alpha))/I)^(1/2));

mf1 =-dDdomegax*(mxa2/Ixd) + dDdalpha*(2*omegaa*ma1/Fa);
% mf1 =-dDdomegax*(mxa2/Ixd) - dDdalpha*(2*omegaa*ma1/Fa);
mf2 = dDdomegax*(mxa1/Ixd) + dDdalpha*(2*omegaa*ma2/Fa);
% mf2 = dDdomegax*(mxa1/Ixd) - dDdalpha*(2*omegaa*ma2/Fa);
b = sqrt(mf1^2 + mf2^2);
diff(b,alpha)







