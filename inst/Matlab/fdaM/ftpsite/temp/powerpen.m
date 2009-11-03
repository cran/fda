function penaltymat = powerpen(basisobj, Lfdobj)
%  POWERPEN  Computes the power basis penalty matrix.
%  Arguments:
%  BASISFD ... a power basis object
%  LFDOBJ  ... either the order of derivative or a
%               linear differential operator to be penalized.
%  Returns a list the first element of which is the basis matrix
%   and the second element of which is the diagonal of the penalty matrix.

%  Last modified:  20 July 2006

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

%  get highest order of derivative and check

nderiv = getnderiv(Lfdobj);
if (nderiv < 0)
    error('NDERIV is negative.');
end

%  compute penalty matrix


if isinteger(Lfdobj)
    
    %  compute exactly if operator is D^NDERIV
    
    exponents = getbasispar(basisobj);
    if length(Lfdobj) == 1
        nderiv = round(Lfdobj);
        if nderiv ~= Lfdobj
            error('Order of derivative must be an integer');
        end
        if nderiv < 0
            error('Order of derivative must be 0 or positive.');
        end
    else
        error('Order of derivative must be a single number.');
    end
    if any(exponents - nderiv < 0) && rang(1) == 0
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
                ifac = ifac*(ideg - k + 1);
            end
        end
        for jbasis=1:ibasis
            jdeg = exponents(jbasis);
            if nderiv == 0
                jfac = 1;
            else
                jfac = jdeg;
                for k=2:nderiv
                    jfac = jfac*(jdeg - k + 1);
                end
            end
            if ideg >= nderiv && jdeg >= nderiv
                penaltymat(ibasis,jbasis) = ifac*jfac* ...
                    (xrange(2)^(ideg+jdeg-2*nderiv+1) -  ...
                    xrange(1)^(ideg+jdeg-2*nderiv+1));
                penaltymat(jbasis,ibasis) = penaltymat(ibasis,jbasis);
            end
        end
    end
else
    
    %  LFDOBJ is not D^NDERIV, use approximate integration by calling
    %  function INPROD().
    
    penaltymat = inprod(basisobj, basisobj, Lfdobj, Lfdobj);
    
end
