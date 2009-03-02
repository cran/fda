function penaltymatrix = polygpen(basisobj, Lfdobj)
%POLYGPEN computes the polygonal penalty matrix for penalty LFD.
%  Arguments:
%  BASISOBJ ... a basis object of type 'bspline'
%  LFDOBJ   ... A linear differential operator object.
%          The highest derivative must be either 0 or 1.
%  Returns the penalty matrix.

%  Last modified:  23 January 2004

%  check BASISOBJ

if ~isa_basis(basisobj)
    error('First argument is not a basis object.');
end

%  check basis type

type = getbasistype(basisobj);
if ~strcmp(type, 'polyg')
    error('BASISOBJ not of type polyg');
end

%  set up default linear differential operator

if nargin < 2, 
    Lfdobj = int2Lfd(0); 
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

%  get highest order of derivative and check

nderiv = getnderiv(Lfdobj);
if (nderiv < 0)
    error('NDERIV is negative.');
end
if (nderiv > 1)
    error('Derivative greater than 1 cannot be taken for polygonal basis.');
end

%  compute penalty matrix

if isinteger(Lfdobj)
    
    %  compute matrix exactly when the operator is D^NDERIV
    
    args    = getbasispar(basisobj);
    n       = length(args);
    argdiff = diff(args);
    penaltymatrix = eye(n);
    if (getnderiv(Lfdobj) == 1)
        penaltymatrix(1,1) = argdiff(  1)/3;
        penaltymatrix(n,n) = argdiff(n-1)/3;
        indx = 2:(n-1);
        for i=indx
            penaltymatrix(i,i) = (argdiff(i)+argdiff(i-1))./3;
        end
        indx = 2:n;
        for i=indx
            penaltymatrix(i,i-1) = argdiff(i-1)/6;
            penaltymatrix(i-1,i) = argdiff(i-1)/6;
        end
    else
        argdiff = 1./argdiff;
        penaltymatrix(1,1) = argdiff(  1);
        penaltymatrix(n,n) = argdiff(n-1);
        indx = 2:(n-1);
        for i=indx
            penaltymatrix(i,i) = argdiff(i)+argdiff(i-1);
        end
        indx = 2:n;
        for i=indx
            penaltymatrix(i,i-1) = -argdiff(i-1);
            penaltymatrix(i-1,i) = -argdiff(i-1);
        end
    end
else
    
    %  LFDOBJ is not D^NDERIV, use approximate integration by calling
    %  function INPROD().
    
    penaltymatrix = inprod(basisobj, basisobj, Lfdobj, Lfdobj);
    
end

penaltymatrix = sparse(penaltymatrix);

