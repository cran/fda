function [res, Dres] = CSTRfitLS(coef, datstruct, fitstruct, ...
                                 lambda, gradwrd)

%  Last modified 9 May 2005

if nargin < 5, gradwrd = 0;  end

fit = fitstruct.fit;

[N, nbasis] = size(datstruct.basismat);
ind1    = 1:nbasis;
ind2    = (nbasis+1):(2*nbasis);
onesb   = ones(1,nbasis);
zeromat = zeros(N, nbasis);

Ccoef  = coef(ind1);
Tcoef  = coef(ind2);

Cwt  = datstruct.Cwt;
Twt  = datstruct.Twt;

Sres = [];
yobs = datstruct.y;
basismat = datstruct.basismat;

if fit(1)
    resC = yobs(:,1) - basismat*Ccoef;
    Sres = [Sres; resC./sqrt(Cwt)];
end

if fit(2)
    resT = yobs(:,2) - basismat*Tcoef;
    Sres = [Sres; resT./sqrt(Twt)];
end

kref   = fitstruct.kref;
EoverR = fitstruct.EoverR;
a      = fitstruct.a;
b      = fitstruct.b;

%  basis function values at quadrature points

quadmat  = datstruct.quadbasismat;
Dquadmat = datstruct.Dquadbasismat;
[nquad, nbasis] = size(quadmat);
onesb = ones(1,nbasis);
onesq = ones(nquad, 1);

%  set up the values of C and T at quad. pts.

Chatquad = quadmat*Ccoef;
Thatquad = quadmat*Tcoef;

DC = Dquadmat*Ccoef;
DT = Dquadmat*Tcoef;

%  set up some constants that are required

V      = fitstruct.V;
rho    = fitstruct.rho;
rhoc   = fitstruct.rhoc;
delH   = fitstruct.delH;
Cp     = fitstruct.Cp;
Cpc    = fitstruct.Cpc;
Tref   = fitstruct.Tref;

%  these constants can vary.
%  see function CSTR2in for other conditions

Fc  = datstruct.Fc;
F   = datstruct.F;
CA0 = datstruct.CA0;
T0  = datstruct.T0;
Tc  = datstruct.Tcin;

%  compute multipliers of outputs

Tdif    = 1./Thatquad - 1./Tref;
betaCC  = kref.*exp(-1e4.*EoverR.*Tdif);
TCfac   = -delH./(rho.*Cp);
betaTC  = TCfac.*betaCC;
aFc2b   = a.*Fc.^b;
K1      = V.*rho.*Cp;
K2      = 1./(2.*rhoc.*Cpc);
betaTT  = Fc.*aFc2b./(K1.*(Fc + K2.*aFc2b));
betaTT0 = F./V;

%  compute right sides of equations 

DChat = -(betaTT0 + betaCC).*Chatquad + betaTT0.*CA0; 
DThat = -(betaTT0 + betaTT).*Thatquad + betaTC.*Chatquad ...
              + betaTT0.*T0 + betaTT.*Tc;

LC = DC - DChat;
LT = DT - DThat;

quadwts = datstruct.quadwts;
rootwts = sqrt(quadwts);

lambdaC = lambda(1);
lambdaT = lambda(2);

Lres = [LC.*rootwts.*sqrt(lambdaC./Cwt); ...
        LT.*rootwts.*sqrt(lambdaT./Twt)];

res  = [Sres; Lres];

%  compute gradient if required

if gradwrd

    %  Derivatives of fit residuals
    
    DSres = [];
    if fit(1)
        DSres = [DSres; [-basismat./sqrt(Cwt), zeromat]];
    end
    if fit(2)
        DSres = [DSres; [zeromat, -basismat./sqrt(Twt)]];
    end

    %  Derivatives of weight functions
    
    DtbetaCC = (1e4.*EoverR./Thatquad.^2).*betaCC;
    DtbetaTC = TCfac.*DtbetaCC;
    
    %  Derivatives of RHS of operators
    
    DcDChat  = -(betaCC + betaTT0);
    DtDChat  = -DtbetaCC.*Chatquad;
    DcDThat  =  betaTC;
    DtDThat  = -(betaTT+betaTT0) + DtbetaTC.*Chatquad;
    
    %  Operator derivatives
    
    DcLC   = Dquadmat - (DcDChat*onesb).*quadmat;
    DtLC   =          - (DtDChat*onesb).*quadmat;
    DcLT   =          - (DcDThat*onesb).*quadmat;
    DtLT   = Dquadmat - (DtDThat*onesb).*quadmat;
    
    %  Multiply operator derivatives by root of
    %  quadrature weights over root of SSE weights
    
    wtmat   = rootwts*onesb;
    DcLCmat = DcLC.*wtmat.*sqrt(lambdaC./Cwt);
    DtLCmat = DtLC.*wtmat.*sqrt(lambdaC./Cwt);
    DcLTmat = DcLT.*wtmat.*sqrt(lambdaT./Twt);
    DtLTmat = DtLT.*wtmat.*sqrt(lambdaT./Twt);
    
    %  Matrices of derivative of operator residuals
    
    DLres  = [[DcLCmat, DtLCmat]; ...
              [DcLTmat, DtLTmat]];
          
    %  Combine with derivative of fit residuals
    
    Dres   = [DSres; DLres];
    
else
    
    Dres = [];
    
end
