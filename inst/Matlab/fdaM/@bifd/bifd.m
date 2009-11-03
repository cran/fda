function bifdobj = bifd(coef, sbasisobj, tbasisobj, bifdnames)
%  BIFD  Creates a bi-functional data object.
%  A bi-functional datum
%  object consists of two bases for expanding a bivariate function and
%  a set of coefficients defining this expansion.  Each basis is contained
%  in a 'basis' object.  That is, a realization of the 'basis' class.
%  Arguments
%  COEF     ... a two-, three-, or four-dimensional array containing
%               coefficient values for the expansion of each set of bivariate
%               function values=terms of a set of basis function values
%               If COEF is a two-way, it is assumed that there is only
%                 one variable and only one replication, and then
%                 the first and second dimensions correspond to
%                 the basis functions for the first and second argument,
%                 respectively.
%               If COEF is a three-way, it is assumed that there is only
%                 one variable per replication, and then
%                 the first and second dimensions correspond to
%                 the basis functions for the first and second argument,
%                 respectively, and the third dimension corresponds to
%                 replications.
%               If COEF is a four-way array, then the fourth dimension
%                 corresponds to variables
%  SBASISOBJ ... a functional data basis object for the first  argument s
%  TBASISOBJ ... a functional data basis object for the second argument t
%  BIFDNAMES ... A cell of length 3 with members containing
%               1. a single name for the argument domain, such as 'Time'
%               2. a name for the replications or cases
%               3. a name for the function.
%  Returns
%  BIFDOBJ  ... a bivariate functional data object

%  last modified 20 July 2006

  superiorto('double', 'struct', 'cell', 'char', 'inline', ...
             'basis');

  if nargin == 0
    bifdobj.coef      = zeros(1,1);
    bifdobj.sbasisobj = basis('bspline',[0,1],3,0.5);
    bifdobj.tbasisobj = basis('bspline',[0,1],3,0.5);
    bifdnames{1} = 'time';
    bifdnames{2} = 'reps';
    bifdnames{3} = 'values';
    bifdobj.fdnames  = bifdnames;
    bifdobj = class(bifdobj, 'bifd');
    return;
  end

  if isa(coef, 'bifd')
    bifdobj = coef;
    return; 
  end

  if nargin < 4
    bifdnames{1} = 'time';
    bifdnames{2} = 'reps';
    bifdnames{3} = 'values';
  end

  coefd = size(coef);
  ndim  = length(coefd);
  if ndim < 2 || ndim > 4
    error(' First argument not of dimension 2, 3 or 4');
  end

  if ~isa_basis(sbasisobj)
    error('Argument SBASISOBJ must be of basis class');
  end
  if ~isa_basis(tbasisobj)
    error('Argument TBASISOBJ must be of basis class');
  end
  if coefd(1) ~= getnbasis(sbasisobj)
    error(['First dimension does not match number of basis functions', ...
           'for first argument.']);
  end
  if coefd(2) ~= getnbasis(tbasisobj)
    error(['Second dimension does not match number of basis functions', ...
           'for second argument.']);
  end

  bifdstr.coef      = coef;
  bifdstr.sbasisobj = sbasisobj;
  bifdstr.tbasisobj = tbasisobj;
  bifdstr.bifdnames = bifdnames;

  bifdobj = class(bifdstr, 'bifd');


