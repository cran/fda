function penaltymat = exponpen(basisobj, Lfdobj)
%EXPONPEN computes the exponential penalty matrix for penalty LFD.
%  Arguments:
%  BASISOBJ ... a basis.fd object of type 'expon'
%  LFDOBJ   ... a linear differential operator to be penalized. 
%               Default is 2.
%  Returns PENALTYMAT, the penalty matrix

%  Last modified:  9 September 2003

%  check BASISOBJ

if ~strcmp(class(basisobj), 'basis')
    error('First argument is not a basis object.')
end

%  check basis type

type = getbasistype(basisobj);
if ~strcmp(type, 'expon')
    error ('Basis type is not exponential.');
end

%  set up default linear differential operator

if nargin < 2, 
    Lfdobj = int2Lfd(2); 
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

%  get highest order of derivative and check

nderiv = getnderiv(Lfdobj);
if nderiv < 0, error('NDERIV is negative.');  end

%  compute penalty matrix

if isinteger(Lfdobj)

    %  Compute penalty matrix exactly if the operator is D^NDERIV

    ratevec    = getbasispar(basisobj);
    nrate      = length(ratevec);
    penaltymat = zeros(nrate,nrate);
    rangeval   = getbasisrange(basisobj);
    tl = rangeval(1);
    tu = rangeval(2);
    for irate = 1:nrate
        ratei = ratevec(irate);
        for jrate = 1:irate
            ratej = ratevec(jrate);
            ratesum = ratei + ratej;
            if ratesum ~= 0
                penaltymat(irate,jrate) = (ratei*ratej)^nderiv * ...
                    (exp(ratesum*tu) - exp(ratesum*tl)) / ratesum;
            else
                if nderiv == 0
                    penaltymat(irate,jrate) = tu - tl;
                end
            end
            penaltymat(jrate,irate) = penaltymat(irate,jrate);
        end
    end
else
    
    %  LFDOBJ is not D^NDERIV, use approximate integration by calling
    %  function INPROD().
    
    penaltymat = inprod(basisobj, basisobj, Lfdobj, Lfdobj);
end

