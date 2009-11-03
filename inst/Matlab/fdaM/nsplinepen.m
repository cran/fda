function [penaltymat, iter] = ...
    nsplinepen(basisobj, Lfdobj, rng, sparsewrd)
%NPLINEPEN computes the N-spline (natural spline) penalty matrix for
%penalty LFDOBJ. 
%  Arguments:
%  BASISOBJ  ... a basis object
%  LFDOBJ    ... a linear differential operator object.  
%                Default int2Lfd(2)
%  RNG       ... A range over which the product is evaluated
%  SPARSEWRD ... if 1, return penaltymatrix in sparse storage mode.

%  added by Kris Villez in August 2011 based on bsplinepen file in the
%  FDA toolbox by Jim Ramsay

%  last modified 31 October 2011 by Kris Villez: changed terminology

%  check BASISOBJ

if ~isa_basis(basisobj)
    error('BASISOBJ is not a basis object.');
end

%  check basis type

type = getbasistype(basisobj);
if ~strcmp(type, 'nspline')
    error('basisobj not of type natural spline');
end

%  set up default value for SPARSEWRD

if nargin < 4, sparsewrd = 1; end

%  set up default value for RNG

range = getbasisrange(basisobj);
if nargin < 3, rng = range;  end

%  set up default linear differential operator

if nargin < 2, 
    Lfdobj = int2Lfd(2); 
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);
  
%  get basis information

nbasis  = getnbasis(basisobj);
params  = getbasispar(basisobj);

%  default for ITER
iter = 0;

%  if there are no internal knots, use the monomial penalty

if isempty(params)
    basisobj   = create_monomial_basis(range, nbasis, 0:nbasis-1);
    penaltymat = monompen(basisobj, Lfdobj, rng);
    return;
end

%  normal case:  PARAMS is not empty

breaks    = [range(1),params,range(2)];  %  break points
nbreaks   = length(breaks);
ninterval = nbreaks - 1;    %  number of intervals defined by breaks

%  check break values
if length(breaks) < 2
    error('The length of argument breaks is less than 2.');
end

%  get highest order of derivative and check

nderiv = getnderiv(Lfdobj);
if nderiv < 0, error('NDERIV is negative'); end
norder = nbasis - length(params) +2;

%  check for order of derivative being equal or greater than
%  order of spline
if nderiv >= norder
    error(['Derivative of order ', num2str(nderiv),                    ...
           ' cannot be taken for N-spline of order ', num2str(norder), ...
           ' Probable cause is a value of the NBASIS argument in',     ...
           ' function FD that is too small.']);
end

%  check for order of derivative being equal to order of spline
%  minus one, in which case following code won't work.

if nderiv > 0 && nderiv == norder - 1
    error(['Penalty matrix cannot be evaluated for derivative of order ', ...
        num2str(nderiv), ' for N-splines of order ', num2str(norder)]);
end

%  set iter

iter = 0;

%  special case where LFD is D^NDERIV and NDERIV = NORDER - 1

if isinteger(Lfdobj) && nderiv == norder - 1
    %  special case of nderiv = norder - 1
    halfseq    = (breaks(2:nbreaks) + breaks(1:(nbreaks-1)))./2;
    halfmat    = nsplineM(halfseq, breaks, norder, nderiv);
    brwidth    = diff(breaks);
    penaltymat = sparse(halfmat' * diag(brwidth) * halfmat);
    return;
end

%  look for knot multiplicities within the range

intbreaks    = getbasispar(basisobj);
index        = intbreaks >= rng(1) & intbreaks <= rng(2);
intbreaks    = intbreaks(index);
if length(index) > 1
    uniquebreaks = min(diff(intbreaks)) > 0;
else
    uniquebreaks = 1;
end

%  if LFD is D^NDERIV, and there are no break multiplicities,
%  use exact computation
if isinteger(Lfdobj) && rng(1) == range(1) && ...
                        rng(2) == range(2) && uniquebreaks
  
    % changes for N-splines:
    % Matrix is computed as for B-splines. Projection is used to convert to
    % N-spline result at the end
                    
    %  Set up the knot sequence

    knots = [range(1)*ones(1,norder), ...
             breaks(2:(nbreaks-1)),          ...
             range(2)*ones(1,norder)]; 

    % Construct  the piecewise polynomial representation
    
    polyorder = norder - nderiv                     ;
    ndegree   = polyorder - 1                       ;
    prodorder = 2*ndegree + 1                       ;   % order of product
    polycoef  = zeros(ninterval, polyorder, norder) ; 
    indxdown  = norder:-1:nderiv+1                  ;
    for i = 1:nbasis +2
        %  compute polynomial representation of B(i,norder,t)(x)
        Coeff    = ppBspline(knots(i:i+norder));
        nrowcoef = size(Coeff,1);
        onescoef = ones(nrowcoef,1);
        % convert the index of the breaks in t to the index in the
        % variable 'breaks'
        CoeffD = Coeff(:,1:polyorder);
        if nderiv > 0
            for ideriv=1:nderiv
                fac = indxdown - ideriv;
                CoeffD = (onescoef*fac).*CoeffD;
            end
        end
        % add the polynomial representation of B(i,norder,t)(x) to f
        if i >= norder, k = norder;  else k = i;       end
        if i <= norder, m = i;       else m = norder;  end
        for j=1:nrowcoef
            polycoef(i-k+j,:,m-j+1) = CoeffD(j,:);
        end
    end
    
    % Compute the scalar products 
    
    prodmat = zeros(nbasis+2);
    convmat = zeros(norder, norder, prodorder);
    for in = 1:ninterval
        %  get the coefficients for the polynomials for this interval
        Coeff = squeeze(polycoef(in,:,:));
        %  compute the coefficients for the products
        for i=0:ndegree-1
            ind = (0:i) + 1;
            convmat(:,:,i+1        ) = ...
                Coeff(ind,          :)'*Coeff(i-ind+2,      :);
            convmat(:,:,prodorder-i) = ...
                Coeff(ndegree-ind+2,:)'*Coeff(ndegree-i+ind,:);
        end
        ind = (0:ndegree)+1;
        convmat(:,:,ndegree+1) = Coeff(ind,:)'*Coeff(ndegree-ind+2,:);
        %  compute the coefficients of the integral
        delta    = breaks(in+1) - breaks(in);
        power    = delta;
        prodmati = zeros(norder);
        for i=1:prodorder
            prodmati = prodmati + ...
                power.*squeeze(convmat(:,:,prodorder-i+1))./i;
            power = power*delta;
        end
        % add the integral to s
        index = in:in+norder-1;
        prodmat(index,index) = prodmat(index,index) + prodmati; 
    end
    
    % projection
    [~,P]       =   nsplineM(breaks, breaks, norder, nderiv, sparsewrd);
    prodmat     =   P'*prodmat*P ;
    
else
    
    % changes for N-splines:
    %   None: inprod does the trick
                    
    %  LFDOBJ is not D^NDERIV, use approximate integration by calling
    %  function INPROD().
    if uniquebreaks
        %  no knot multiplicities
        [prodmat, iter] = inprod(basisobj, basisobj, Lfdobj, Lfdobj, rng);
    else
        %  knot multiplicities ... find their locations
        rngvec = rng(1);
        for i=2:nbreaks
            if breaks(i) == breaks(i-1)
                rngvec = [rngvec, breaks(i)];
            end
        end
        rngvec = unique(rngvec);
        nrng   = length(rngvec);
        if rngvec(nrng) < rng(2)
            rngvec = [rngvec,rng(2)];
            nrng   = nrng + 1;
        end
        %  sum prodmat over intervals between knot multiplicities
        prodmat = zeros(nbasis);
        for i = 2:nrng
            rngi = [rngvec(i-1) + 1e-10, rngvec(i) - 1e-10];
            [prodmati, iter] = ...
                 inprod(basisobj, basisobj, Lfdobj, Lfdobj, rngi);
            prodmat = prodmat + prodmati;
        end
    end
    
end

if sparsewrd
    penaltymat = sparse(prodmat);
else
    penaltymat = prodmat;
end
