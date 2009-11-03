function basisobj = create_TP_basis(basisobj1, basisobj2)
%  CREATE_TP_BASIS sets up a tensor product basis for the analysis
%  of two-dimensional spatial data or of time crossed with one-dimensional 
%  data.  It constructs this basis from two one-dimensional basis objects.
%
%  Arguments:
%  BASISOBJ1 ... A one-dimensional functional basis object.
%  BASISOBJ2 ... A one-dimensional functional basis object.
%  
%  Returns:
%  An object of the basis class with parameters stored in member PARAMS.
%  In this case the member PARAMS is a struct object with members
%  BASISOBJ1 and BASISOBJ2.  The RANGEVAL object is a 2 by 2 matrix with
%  the first row containing the lower and upper range values for object
%  BASISOBJ1 and the second row the same for BASISOBJ2.   Similarly,
%  the NBASIS member is a vector of length 2 containing the two numbers
%  of basis functions.   The TYPE member is 'TP'.

%  Functional data objects are defined

%  Last modified 7 March 2011 by Jim Ramsay

%  check for number of arguments

if nargin < 2
    error('Less than two input arguments.');
end

%  check the two arguments

if ~strcmp(class(basisobj1), 'basis') || ...
   ~strcmp(class(basisobj2), 'basis')
    error('At least one of the arguments is not of the basis class.')
end

%  set the TYPE member

type = 'TP';

%  member RANGEVAL is empty, and access to the two ranges is obtained
%  through the PARAMS member which contains the two one-dimensional
%  basis objects.

rangeval = [];

% member NBASIS is the product of the two numbers of basis functions.
% evaluation of a functional data object defined on a basis of the TP
% class will expect a vector of coefficients of this size, and will
% convert these to a matrix of size NBASIS1 by NBASIS2 by the reshape
% function.

nbasis = getnbasis(basisobj1)*getnbasis(basisobj2);

%  member PARAMS is a struct object containing the two basis objects
%  as members

params.basis1 = basisobj1;
params.basis2 = basisobj2;

%  define the TP basis object

basisobj = basis(type, rangeval, nbasis, params);
