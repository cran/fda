function [Dy, DDy] = CSTR2(tobs, y, fitstruct, condition, Tlim)
%  CSTR2 is called by the ODE solving function and evaluates the
%  right hand side of the equation.

%  Last modified 5 May 2005

tau = 1;

[F,CA0,T0,Tcin,Fc] = CSTR2in(tobs, condition, tau);

kref   = fitstruct.kref;
EoverR = fitstruct.EoverR;
a      = fitstruct.a;
b      = fitstruct.b;

V      = fitstruct.V;
FoverV = F./V;

aFc2bp = a.*Fc.^b;
const1 = V.*fitstruct.rho.*fitstruct.Cp;
const2 = 2.*fitstruct.rhoc.*fitstruct.Cpc;

betaTTcin = Fc.*aFc2bp./(const1.*(Fc+aFc2bp./const2));
betaTT0   = FoverV;

%  set up current values of variables

if length(y) == 2 
    C = y(1);
    T = y(2);
else
    C = y(:,1);
    T = y(:,2);
end

Tdif = 1./T - 1/fitstruct.Tref;
betaCC = kref*exp(-1e4*EoverR*Tdif);
betaTC = (-fitstruct.delH/const1).*betaCC;
betaTT = FoverV + betaTTcin;

DC = -(FoverV + betaCC).*C  + (F./V).*CA0;
DT = -betaTT.*T + betaTC.*C + (F./V).*T0 + betaTTcin.*Tcin;

if length(y) == 2
    Dy  = [DC; DT];
else
    Dy  = [DC, DT];
end

waitbar(tobs/Tlim)