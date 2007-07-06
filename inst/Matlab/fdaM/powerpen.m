function penaltymat = powerpen(basisobj, Lfdobj)
%  POWERPEN  Computes the power basis penalty matrix.
%  Arguments:
%  BASISFD ... a power basis object
%  LFDOBJ  ... either the order of derivative or a
%               linear differential operator to be penalized.

%  Last modified:  4 January 2008

%  check BASISOBJ

if ~strcmp(class(basisobj), 'basis')
    error('First argument is not a basis object.');
end

%  check basis type

type = getbasistype(basisobj);
rang = getbasisrange(basisobj);
if ~strcmp(type, 'power')
    error('BASISOBJ not of type POWER.');
end

%  set up default linear differential operator

if nargin < 2, 
    Lfdobj = int2Lfd(2); 
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

%  compute penalty matrix

if isinteger(Lfdobj)
    
    %  compute exactly if operator is D^NDERIV
    
    nderiv = getnderiv(Lfdobj);
    exponents = getbasispar(basisobj);
    if any(exponents - nderiv < 0) && rang(1) <= 0
        error('A negative exponent is needed and an argument value is 0.');
    end
    nbasis = getnbasis(basisobj);
    penaltymat = zeros(nbasis);
    xrange = getbasisrange(basisobj);
    for ibasis=1:nbasis
        ideg = exponents(ibasis);
        if nderiv == 0
            ifac = 1;
        else
            ifac = ideg;
            for k=2:nderiv
                if ideg == k-1
                    ifac = 0;
                else
                    ifac = ifac*(ideg - k + 1);
                end
            end
        end
        for jbasis=1:ibasis-1
            jdeg = exponents(jbasis);
            if nderiv == 0
                jfac = 1;
            else
                jfac = jdeg;
                for k=2:nderiv
                    if jdeg == k-1
                        jfac = 0;
                    else
                        jfac = jfac*(jdeg - k + 1);
                    end
                end
            end
            if ifac*jfac == 0
                penaltymat(ibasis,jbasis) = 0;
            else
                penaltymat(ibasis,jbasis) = ifac*jfac* ...
                    (xrange(2)^(ideg+jdeg-2*nderiv+1) -  ...
                     xrange(1)^(ideg+jdeg-2*nderiv+1))/ ...
                    (ideg + jdeg - 2*nderiv + 1);
            end
            penaltymat(jbasis,ibasis) = penaltymat(ibasis,jbasis);
        end
        if ifac == 0
            penaltymat(ibasis,ibasis) = 0;
        elseif 2*ideg - 2*nderiv + 1 == 0
            penaltymat(ibasis,ibasis) = ifac^2* ...
                    (log(xrange(2)) - log(xrange(1)));
        else
            penaltymat(ibasis,ibasis) = ifac^2* ...
                    (xrange(2)^(2*ideg-2*nderiv+1) -  ...
                     xrange(1)^(2*ideg-2*nderiv+1))/ ...
                     (2*ideg - 2*nderiv + 1);
        end
    end
else
    
    %  LFDOBJ is not D^NDERIV, use approximate integration by calling
    %  function INPROD().
    
    penaltymat = inprod(basisobj, basisobj, Lfdobj, Lfdobj);
    
end

