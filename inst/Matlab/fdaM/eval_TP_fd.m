function surfarray = eval_TP_fd(Xvec, Yvec, fdobj, xderiv, yderiv)
% EVAL_TP_FD evaluates the TP fd object at points (Xvec,Yvec)
%
%        arguments:
% FELSPLOBJ a FELspline object
% XVEC   ... an array of x-coordinates.
% YVEC   ... an array of y-coordinates.
% XDERIV ... order of derivative with respect to x
% YDERIV ... order of derivative with respect to y
% 
%        output:
% SURFARRAY an array of the same size as Xvec and Yvec containing the value  
%           of FDOBJ at (Xvec,Yvec).

%  Last modified on 7 March 2011 by Jim ramsay.

%  Set up the arguments if the first argument is a matrix with two
%  columns

if size(Xvec,2) == 2
    if nargin < 2
        error(['First argument is a coordinate matrix and ', ...
               'the second argument is not supplied.']);
    end
    if nargin < 4,  xderiv = 0;  end
    if nargin < 3,  yderiv = 0;  end
    fdobj = Yvec;
    Yvec  = Xvec(:,2);
    Xvec  = Xvec(:,1);
else
    if nargin < 3
        error(['First and second arguments are coordinate vectors ', ...
               'and the third argument is not supplied.']);
    end
    if nargin < 5,  xderiv = 0;  end
    if nargin < 4,  yderiv = 0;  end
end

%  check Xvec

if ~isa(Xvec,'double')
   error('Xvec is not a numerical array')
else
   Xvec  = Xvec(:);     % treat Xvec as a column vector
end

%  check Yvec

if ~isa(Yvec,'double')
   error('Yvec is not a numerical array')
else
   Yvec=Yvec(:);     % treat Yvec as a column vector
end

N1 = length(Xvec);
N2 = length(Yvec);

%  check the type of FDOBJ

basisobj = getbasis(fdobj);
if ~strcmp(getbasistype(basisobj), 'TP')
    error('The basis object for FDOBJ is not of type TP.');
end

%  extract the two basis objects

basis_struct = getbasispar(basisobj);
basisobj1 = basis_struct.basis1;
basisobj2 = basis_struct.basis2;

%  extract the two numbers of basis functions

nbasis1 = getnbasis(basisobj1);
nbasis2 = getnbasis(basisobj2);

%  extract the two ranges

rangeval1 = getbasisrange(basisobj1);
rangeval2 = getbasisrange(basisobj2);

%  evaluate the two basis objects:

basismat1 = eval_basis(Xvec, basisobj1, xderiv);
basismat2 = eval_basis(Yvec, basisobj2, yderiv);

%  extract the coefficient array

coef = getcoef(fdobj);

coefdim = size(coef);
nbasis = coefdim(1);
if nbasis ~= nbasis1*nbasis2
    error(['First dimension of coefficient array not equal to ', ...
           'product of two numbers of basis functions.']);
end
ndim  = length(coefdim);
nsurf = coefdim(2);
if ndim == 2
    nvar = 1;
else
    nvar = coefdim(3);
end

%  loop through variables and surfaces

if ndim == 2
    surfarray = zeros(nbasis1,nbasis2,nsurf);
else
    surfarray = zeros(nbasis1,nbasis2,nsurf,nvar);
end

for ivar = 1:nvar
    for isurf = 1:nsurf
        if ndim == 2
            coefvec = coef(:,isurf);
            coefmat = reshape(coefvec,nbasis1,nbasis2);
            surfmat = basismat1*coefmat*basismat2';
            if nsurf > 1
                surfarray(:,:,isurf) = surfmat;
            else
                surfarray = surfmat;
            end
        else
            coefvec = coef(:,isurf,ivar);
            coefmat = reshape(coefvec,nbasis1,nbasis2);
            surfmat = basismat1*coefmat*basismat2';
            surfarray(:,:,isurf,ivar) = surfmat;
        end            
    end
end
    



