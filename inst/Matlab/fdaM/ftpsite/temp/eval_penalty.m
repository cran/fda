function [penaltymat, iter] = eval_penalty(basisobj, Lfdobj, rng)
%  EVAL_PENALTY evaluates the inner products of a linear 
%  differential operator L defined by LFDOBJ applied to a set of 
%  basis functions defined by BASISOBJ.
%
%  LFDOBJ is a functional data object defining the order m 
%  NONHOMOGENEOUS linear differential operator of the form
%  Lx(t) = w_0(t) x(t) + ... + w_{m-1}(t) D^{m-1}x(t) + 
%          \exp[w_m(t)] D^m x(t) + ...
%          a_1(t) u_1(t)  + ... + a_k(t) u_k(t).
%  This is a change from previous usage where LFDOBJ was assumed to 
%  define a HOMOGONEOUS differential operator.  See function
%  @Lfd/Lfd() for details.
%
%  Arguments:
%  BASISOBJ ... Either a basis object or an fd object or 
%               an fdPar object.  If an fdPar object,
%               and if no LFDOBJ is supplied, the LFDOBJ
%               in the fdPar object is used.
%  LFDOBJ   ... A linear differential operator object
%               applied to the functions that are evaluated.
%  RNG      ... A range over which the product is evaluated
%
%  Returns:  
%  PENALTYMAT ... Symmetric matrix containing inner products.
%                 This matrix should be non-negative definite 
%                 With NDERIV zero eigenvalues, where NDERIV
%                 is the highest order derivative in LFDOBJ.
%                 However, rounding error will likely cause 
%                 NDERIV smallest eigenvalues to be nonzero,
%                 so be careful about calling CHOL or otherwise
%                 assuming the range is N - NDERIV.  
%  ITER       ... Number of iterations taken in INPROD if
%                 it is called.  Otherwise 0.
%                

%  last modified 20 July 2006

%  set ITER to default value of 0

iter = 0;

%  check BASISOBJ

if isa_fd(basisobj)
    basisobj = getbasis(basisobj);
end

Lfd_default = int2Lfd(0);

if isa_fdPar(basisobj)
    if nargin < 2
        %  if the penalty matrix is already stored in the
        %  fdPar object, just get it and return it.
        penaltymat = getpenmat(basisobj);
        if ~isempty(penaltymat)
            return
        end
    end
    Lfd_default = getLfd(basisobj);
    basisobj    = getbasis(getfd(basisobj));
end
if ~isa_basis(basisobj)
    error('Argument BASISOBJ is not a functional basis object.');
end

%  set up default values

if nargin < 3, rng = getbasisrange(basisobj);  end
if nargin < 2, Lfdobj = Lfd_default;           end

%  deal with the case where LFDOBJ is an integer

if isnumeric(Lfdobj)
    nderiv = Lfdobj;
    if nderiv ~= round(Lfdobj)
        error('LFDOBJ numeric but not an integer.');
    end
    if nderiv < 0
        error('LFDOBJ an integer but negative.');
    end
    Lfdobj = int2Lfd(nderiv);
end
    
%  check LFDOBJ

if ~isa_Lfd(Lfdobj)
    error (['Argument LFDOBJ is not a ', ...
            'linear differential operator object.']);
end

%  determine basis type

type = getbasistype(basisobj);

%  choose appropriate penalty matrix function

switch type
    case 'bspline'
        [penaltymat, iter] = bsplinepen(basisobj, Lfdobj, rng);
    case 'const'
        rangeval   = getbasisrange(basisobj);
        penaltymat = rangeval(2) - rangeval(1);
    case 'expon'
        penaltymat = exponpen(basisobj,   Lfdobj);
    case 'fourier'
        penaltymat = fourierpen(basisobj, Lfdobj);
    case 'monom'
        penaltymat = monompen(basisobj,   Lfdobj);
    case 'polyg'
        penaltymat = polygpen(basisobj,   Lfdobj);
    case 'power' 
        penaltymat = powerpen(basisobj,   Lfdobj);   otherwise
        error('Basis type not recognizable');
end

%  If drop indices are provided, drop rows and cols 
%  associated with these indices

dropind = getdropind(basisobj);
nbasis  = getnbasis(basisobj);
if length(dropind) > 0
    index = 1:nbasis;
    for i=1:length(dropind)
        index = index(index ~= dropind(i));
    end
    penaltymat = penaltymat(index,index);
end

%  Make matrix symmetric since small rounding errors can
%  sometimes results in small asymmetries

penaltymat = (penaltymat + penaltymat')./2;


