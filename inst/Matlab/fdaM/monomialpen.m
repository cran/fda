function penaltymat = monomialpen(basisobj, Lfdobj, rng)
%  MONOMPEN  Computes the monomial penalty matrix.
%
%  Warning:  This version is incomplete in the sense that
%    it only works with LFDOBJ = D^m
%
%  Arguments:
%  BASISFD  ... a monomial basis object
%  LFDOBJ   ... either the order of derivative or a
%               linear differential operator to be penalized.
%  Returns a list the first element of which is the basis matrix
%   and the second element of which is the diagonal of the penalty matrix.

%  Last modified:  1 November 2007

%  check BASISOBJ

if ~strcmp(class(basisobj), 'basis')
    error('First argument is not a basis.fd object.');
end
type = getbasistype(basisobj);
if ~strcmp(type, 'monom')
    error('BASISOBJ not of type monom');
end

%  set default value of LFDOBJ

if nargin < 2
    Lfdobj = int2Lfd(2);
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

%  check whether LFDOBJ is of form D^m

if ~isinteger(Lfdobj)
    error('This version cannot handle noninteger operators.');
end

nderiv = getnderiv(Lfdobj);

%  get basis information

nbasis     = getnbasis(basisobj);
if nargin < 3
    xrange = getbasisrange(basisobj);
else
    xrange = rng;
end
exponents  = getbasispar(basisobj);

%  compute the penalty matrix

penaltymat = zeros(nbasis);
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
            ipow = ideg+jdeg-2*nderiv+1;
            penaltymat(ibasis,jbasis) = ifac*jfac* ...
                (xrange(2)^ipow - xrange(1)^ipow)/ipow;
            penaltymat(jbasis,ibasis) = penaltymat(ibasis,jbasis);
        end
    end
end

%  If drop indices are provided, drop rows and cols
%  associated with these indices

dropind = getdropind(basisobj);
if ~isempty(dropind)
    nbasis  = getnbasis(basisobj) + length(dropind);
    index = 1:nbasis;
    for i=1:length(dropind)
        index = index(index ~= dropind(i));
    end
    penaltymat = penaltymat(index,index);
end


