function [res, Dres, fitstruct, df, gcv] = ...
    CSTRfn(parvec, datstruct, fitstruct, CSTRbasis, lambda, gradwrd)

%  Last modified 29 July 2005
%  Cosmetic modifications 2007.07.21 by Spencer Graves
%  commenting out variables (F1, gradnorm0, rootwts) never used 

if nargin < 6, gradwrd = 1;  end

fit = fitstruct.fit;
%%
% 1.  load parameters

[kref, EoverR, a, b] = par2vals(parvec, fitstruct);

%  set up fitstruct

fitstruct.kref   = kref;
fitstruct.EoverR = EoverR;
fitstruct.a      = a;
fitstruct.b      = b;

estimate = fitstruct.estimate;
%%
% 2.  Inner optimization:  optimize fit with respect to coef

tolval  = 1e-10;
itermax = 10;
coef0   = fitstruct.coef0;
ncoef   = length(coef0);
%  initial value
[res0, Dres] = CSTRfitLS(coef0, datstruct, fitstruct, ...
                         lambda, 1);
F0 = mean(res0.^2) ;
% F1 = 0;
iter = 0;
fundif = 1;
% gradnorm0 = mean(res0'*Dres);
gradnorm1 = 1;
%%
% 3.  Gauss-Newton optimization loop
while ((gradnorm1 > tolval) || (fundif > tolval))
    iter = iter + 1;
    if iter > itermax, break; end
    Dcoef = Dres\res0;
    %  initial step:  alpha = 1
    coef1 = coef0 - Dcoef;
    [res1, Dres] = CSTRfitLS(coef1, datstruct, fitstruct, ...
        lambda, 1);
    F1 = mean(res1.^2);
%     alpha = 1;
%     %  smaller steps as required, halving alpha each time
%     while F1 >= F0
%         alpha = alpha/2;
%         coef1 = coef0 - alpha.*Dcoef;
%         [res1, Dres] = CSTRfitLS(coef1, datstruct, fitstruct, ...
%             lambda, 1);
%         F1 = mean(res1.^2);
%     end
% 
% save CSTRfnTst0 -v6 parvec datstruct fitstruct CSTRbasis ... 
% lambda gradwrd res0 F0 gradnorm0 Dcoef F1 ;
% DresFull = full(Dres);  
% dlmwrite('Dres-Matlab.csv', DresFull, ',') ;  
    gradnorm1 = mean(res1'*Dres);
    fundif = abs(F0-F1)/abs(F0);
%     disp(num2str([iter, F1, fundif, gradnorm1]))
    coef0     = coef1;
    res0      = res1;
    F0        = F1;
%    gradnorm0 = gradnorm1;
end
%%
% 4.  update coef

coef = coef0;
fitstruct.coef0 = coef;
%%
% 5.  compute df and gcv

Zmat = Dres(1:ncoef,:);
Nval = length(res1);
Rfac = Dres(ncoef+1:Nval,:);
Smat = Zmat*inv(Zmat'*Zmat + Rfac'*Rfac)*Zmat';
df   = sum(diag(Smat));
dfe  = ncoef - df;
gcv  = (ncoef/dfe)*sum(res1(1:ncoef).^2)/dfe;
%%
% 6.  compute fits and residuals

[N, nbasis] = size(datstruct.basismat);
ind1  = 1:nbasis;
ind2  = (nbasis+1):(2*nbasis);
Ccoef = coef(ind1);
Tcoef = coef(ind2);

phimat = datstruct.basismat;

Chat = phimat*Ccoef;
That = phimat*Tcoef;

yobs = datstruct.y;

Cwt = datstruct.Cwt;
Twt = datstruct.Twt;

res = [];

if fit(1)
    res = [res; (yobs(:,1) - Chat)./sqrt(Cwt)];
end

if fit(2)
    res = [res; (yobs(:,2) - That)./sqrt(Twt)];
end
%%
% 7.  Derivatives? 

if gradwrd
    
    % 7.1.  set up basis and basis derivatve matrices
    
    quadmat  = datstruct.quadbasismat;
    Dquadmat = datstruct.Dquadbasismat;
    
    [nquad, nbasis] = size(quadmat);  
    onesb = ones(1,nbasis);
    onesq = ones(nquad, 1);
    
    % 7.2.  set up some constants that are required
    
    V      = fitstruct.V;
    rho    = fitstruct.rho;
    rhoc   = fitstruct.rhoc;
    delH   = fitstruct.delH;
    Cp     = fitstruct.Cp;
    Cpc    = fitstruct.Cpc;
    Tref   = fitstruct.Tref;

    % 7.3.  Set up input arrays
    
    F   = datstruct.F;   
    CA0 = datstruct.CA0; 
    T0  = datstruct.T0;
    Tc  = datstruct.Tcin;
    Fc  = datstruct.Fc;  

    % 7.4.  C and T values at fine grid
    
    Chat  = quadmat*Ccoef;
    That  = quadmat*Tcoef;
    DChat = Dquadmat*Ccoef;
    DThat = Dquadmat*Tcoef;

    % 7.5.  betaCC and betaTC depend on kref and Eover R
    Tdif   = 1./That - 1./Tref;
    temp   = exp(-1e4.*EoverR.*Tdif);
    betaCC = kref.*temp; 
    TCfac  = -delH./(rho.*Cp);
    betaTC = TCfac.*betaCC;
    % 7.6.  betaTT depends on a and b
    Fc2b    = Fc.^b;
    aFc2b   = a.*Fc2b;
    K1      = V.*rho.*Cp;
    K2      = 1./(2.*rhoc.*Cpc);
    betaTT  = Fc.*aFc2b./(K1.*(Fc + K2.*aFc2b));    
    betaTT0 = F./V;

    % 7.7.  compute derivatives of residuals
    
    %  L values
    
    LC = DChat + (betaTT0 + betaCC).*Chat - betaTT0.*CA0.*onesq;
    LT = DThat + (betaTT0 + betaTT).*That - betaTC.*Chat - ...
                 (betaTT0.*T0 + betaTT.*Tc); 
          
    %  first order derivatives of L values 
    
    %  derivatives of L values with respect to 
    %  coefficient vectors c and t
    
    DtbetaCC = (1e4*EoverR./That.^2).*betaCC;
    DtbetaTC = TCfac.*DtbetaCC;
    
    DcDChat  = -(betaCC + betaTT0);
    DtDChat  = -DtbetaCC.*Chat;
    DcDThat  =  betaTC;
    DtDThat  = -(betaTT+betaTT0) + DtbetaTC.*Chat;
    
    DcLC   = Dquadmat - (DcDChat*onesb).*quadmat;
    DtLC   =          - (DtDChat*onesb).*quadmat;
    DcLT   =          - (DcDThat*onesb).*quadmat;
    DtLT   = Dquadmat - (DtDThat*onesb).*quadmat;

    quadwts    = datstruct.quadwts;
%    rootwts    = sqrt(quadwts);
    quadwtsmat = quadwts*onesb;

    %  k derivatives
    
    lamC = lambda(1);
    lamT = lambda(2);
    
    % 7.8.  assemble the Jacobian matrix
    
    DLC = sqrt(lamC/Cwt).*[DcLC, DtLC];
    DLT = sqrt(lamT/Twt).*[DcLT, DtLT];
    Jacobian = [DLC; DLT];
    
    % 7.9.  compute derivatives with respect to parameters
    
    %  set up right hand side of equation D2GDc
    
    D2GDc = [];
    
    %  kref
    
    if estimate(1)

        %  first derivative of L values

        DkbetaCC = temp;
        DkbetaTC = TCfac.*DkbetaCC;

        DkLC =  DkbetaCC.*Chat;
        DkLT = -DkbetaTC.*Chat;
        
        %  second derivative of L values

        DktbetaCC = (1e4.*EoverR./That.^2).*temp;
        DktbetaTC = TCfac.*DktbetaCC;

        DkcLC = (  DkbetaCC        *onesb).*quadmat;
        DkcLT = ( -DkbetaTC        *onesb).*quadmat;
        DktLC = (( DktbetaCC.*Chat)*onesb).*quadmat;
        DktLT = ((-DktbetaTC.*Chat)*onesb).*quadmat;

        D2GDck = zeros(2*nbasis,1);
        D2GDck(ind1,1) = (lamC/Cwt).* ...
            (DcLC'*(DkLC.*quadwts) + DkcLC'*(LC.*quadwts)) + ...
            (lamT/Twt).* ...
            (DcLT'*(DkLT.*quadwts) + DkcLT'*(LT.*quadwts));
        D2GDck(ind2,1) = (lamC/Cwt).* ...
            (DtLC'*(DkLC.*quadwts) + DktLC'*(LC.*quadwts)) + ...
            (lamT/Twt).* ...
            (DtLT'*(DkLT.*quadwts) + DktLT'*(LT.*quadwts));
        
        D2GDc = [D2GDc, D2GDck];

    end
    
    %  EoverR

    if estimate(2)

        %  first derivative of L values

        Dtemp    = -1e4.*kref.*Tdif.*temp;
        DEbetaCC = Dtemp;
        DEbetaTC = TCfac.*DEbetaCC;

        DELC  =  DEbetaCC.*Chat;
        DELT  = -DEbetaTC.*Chat;

        DEtbetaCC = (1e4.*kref  ./That.^2).* ...
            (1 - 1e4.*EoverR.*Tdif).*temp;
        DEtbetaTC = TCfac.*DEtbetaCC;

        DEcLC = (  DEbetaCC        *onesb).*quadmat;
        DEcLT = ( -DEbetaTC        *onesb).*quadmat;
        DEtLC = (( DEtbetaCC.*Chat)*onesb).*quadmat;
        DEtLT = ((-DEtbetaTC.*Chat)*onesb).*quadmat;
        
        D2GDcE = zeros(2*nbasis,1);
        D2GDcE(ind1,1) = (lamC/Cwt).* ...
            (DcLC'*(DELC.*quadwts) + DEcLC'*(LC.*quadwts)) + ...
            (lamT/Twt).* ...
            (DcLT'*(DELT.*quadwts) + DEcLT'*(LT.*quadwts));
        D2GDcE(ind2,1) = (lamC./Cwt).* ...
            (DtLC'*(DELC.*quadwts) + DEtLC'*(LC.*quadwts)) + ...
            (lamT./Twt).* ...
            (DtLT'*(DELT.*quadwts) + DEtLT'*(LT.*quadwts));
        D2GDc = [D2GDc, D2GDcE];

    end
    
    %  a

    if estimate(3)

        %  first derivative of L values

        DhbetaTT = (betaTT./aFc2b).*(1 - K1.*K2.*betaTT./Fc);
        DabetaTT = DhbetaTT.*Fc2b;
        
        DaLT = DabetaTT.*(That - Tc);
        
        DatLT = (DabetaTT*onesb).*quadmat;
        
        D2GDca = zeros(2*nbasis,1);
        D2GDca(ind1,1) = (lamT/Twt).*(DcLT'*(DaLT.*quadwts));
        D2GDca(ind2,1) = (lamT./Twt).* ...
            (DtLT'*(DaLT.*quadwts) + DatLT'*(LT.*quadwts));
        
        D2GDc = [D2GDc, D2GDca];
        
    end

    if estimate(4)

        %  b derivative of L values

        DhbetaTT = (betaTT./aFc2b).*(1 - K1.*K2.*betaTT./Fc);
        DbbetaTT = DhbetaTT.*b.*aFc2b./Fc;
        
        DbLT = DbbetaTT.*(That - Tc);
        
        DbtLT = (DbbetaTT*onesb).*quadmat;

        D2GDcb = zeros(2*nbasis,1);
        D2GDcb(ind1,1) = (lamT/Twt).*(DcLT'*(DbLT.*quadwts));
        D2GDcb(ind2,1) = (lamT./Twt).* ...
            (DtLT'*(DbLT.*quadwts) + DbtLT'*(LT.*quadwts));
        
        D2GDc = [D2GDc, D2GDcb];
        
    end
%%
% 8.  Construct D2GDc2 
%
% 8.1.  First part 
    Wmat = [quadwtsmat, quadwtsmat; quadwtsmat, quadwtsmat];
    
    D2GDc2  = (Jacobian.*Wmat)'*Jacobian;
    ZtZmat = phimat'*phimat;
    if fit(1)
        D2GDc2(ind1,ind1) = D2GDc2(ind1,ind1) + ZtZmat./Cwt;
    end
    if fit(2)
        D2GDc2(ind2,ind2) = D2GDc2(ind2,ind2) + ZtZmat./Twt;
    end

% 8.2.  Add second derivative information
    
    DttbetaCC = (1e4.*kref.*EoverR./That.^2).* ...
                (1e4.*EoverR./That.^2 - 2./That).*temp;
    DttbetaTC = TCfac.*DttbetaCC;
    
    DctLC = sparse(zeros(nbasis,nbasis));
    DttLC = sparse(zeros(nbasis,nbasis));
    DctLT = sparse(zeros(nbasis,nbasis));
    DttLT = sparse(zeros(nbasis,nbasis));
    
    norder = nbasis - length(getbasispar(CSTRbasis));
    for i=1:nbasis
        jstart = max([1,i-norder+1]);
        for j=jstart:i
            qijvec = quadmat(:,i).*quadmat(:,j).*quadwts;
            DctLC(i,j) = sum(qijvec.*LC.*DtbetaCC);
            if i ~= j, DctLC(j,i) = DctLC(i,j); end
            DttLC(i,j) = sum(qijvec.*LC.*DttbetaCC.*Chat);
            if i ~= j, DttLC(j,i) = DttLC(i,j); end
            DctLT(i,j) = sum(qijvec.*LT.*DtbetaTC);
            if i ~= j, DctLT(j,i) = DctLT(i,j); end
            DttLT(i,j) = sum(qijvec.*LT.*DttbetaTC.*Chat);
            if i ~= j, DttLT(j,i) = DttLT(i,j); end
        end
    end
        
    DctL = lamC.*DctLC./Cwt + lamT.*DctLT./Twt;
    DttL = lamC.*DttLC./Cwt + lamT.*DttLT./Twt;
          
% 8.3.  modify D2GDc2
   
    D2GDc2(ind1,ind2) = D2GDc2(ind1,ind2) + DctL;
    D2GDc2(ind2,ind1) = D2GDc2(ind2,ind1) + DctL';
    D2GDc2(ind2,ind2) = D2GDc2(ind2,ind2) + DttL;
    
% 8.4.  compute (D2GDc2)^{-1} D2GDc
    
    DcDtheta = D2GDc2\D2GDc;
    
% 8.5.  set up Dres
    
    Dres = [];
    if fit(1)
        Dres = [Dres; phimat*DcDtheta(ind1,:)./sqrt(Cwt)];
    end
    if fit(2)
        Dres = [Dres; phimat*DcDtheta(ind2,:)./sqrt(Twt)];
    end
% save CSTRfnD2 res Dres; 
end

%  -------------------------------------------------------

function [kref, EoverR, a, b] = par2vals(parvec, fitstruct)

estimate = fitstruct.estimate;

m = 0;
if estimate(1)
    m = m + 1;
    kref = parvec(m);
else
    kref = fitstruct.kref;
end

if estimate(2)
    m = m + 1;
    EoverR = parvec(m);
else
    EoverR = fitstruct.EoverR;
end

if estimate(3)
    m = m + 1;
    a = parvec(m);
else
    a = fitstruct.a; 
end

if estimate(4)
    m = m + 1;
    b = parvec(m);
else
    b = fitstruct.b; 
end


