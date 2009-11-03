function newbasisobj = putbasisvalues(basisobj, basisvalues)
%  PUTVALUES   Enters a cell array of values of basis functions
%  and a number of their derivatives into basis object BASISOBJ 
%  BASISVALUES{1} contains the argument values associated with
%  the basis values
%  BASISVALUES{2} contains in some form the values of the basis
%  functions, and BASISVALUES{2+IDERIV} contains the corresponding
%  values of the derivative of order IDERIV.
%  For the classic basis systems such as spline, fourier, monomial
%  and etc. BASISVALUES{2+IDERIV}, IDERIV=0,...,NDERIV will be
%  a matrix with number of rows equal to the length of BASISVALUES{1}
%  and number of columns equal to the number of basis functions.
%  However, for some basis systems, this may not apply, and a 
%  specific example is the 'fdVariance' system where the basis
%  value information is in a cell array of length NSURF, the number
%  of surfaces in a compound covariance surface.
%
%  Last modified 18 June 2007

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

%  check values

if ~isempty(basisvalues)
    %  check BASISVALUES
    if ~iscell(basisvalues)
        error('BASISVALUES is not a cell object.')
    end
    sizevec = size(basisvalues);
    if length(sizevec) > 2
        error('BASISVALUES is not 2-dimensional.')
    end
    if strcmp(class(basisvalues{1,2}),'double')
        for i=1:sizevec(1)
            if length(basisvalues{i,1}) ~= size(basisvalues{i,2},1)
                error(['Number of argument valuesdoes not equal ', ...
                    'number of values.']);
            end
        end
    end
end

newbasisobj.type        = basisobj.type;
newbasisobj.rangeval    = basisobj.rangeval;
newbasisobj.nbasis      = basisobj.nbasis;
newbasisobj.params      = basisobj.params;
newbasisobj.dropind     = basisobj.dropind;
newbasisobj.quadvals    = basisobj.quadvals;
newbasisobj.values      = basisobj.values;
newbasisobj.basisvalues = basisvalues;

newbasisobj = class(newbasisobj, 'basis');