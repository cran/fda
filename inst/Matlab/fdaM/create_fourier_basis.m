function basisobj = create_fourier_basis(rangeval, nbasis, period, ...
                                         dropind)
%CREATE_FOURIER_BASIS  Creates a fourier functional data basis.
%  This function is identical to FOURIER_BASIS.
%  Arguments ...
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS   ... the number of basis functions
%  PERIOD   ... The period.  That is, the basis functions are periodic on
%                 the interval (0,PARAMS) or any translation of it.
%  DROPIND  ... a set of indices in 1:NBASIS of basis functions to drop
%                when basis objects are arguments.  Default is [];
%  Returns
%  BASIS_fd  ... a functional data basis object of type 'fourier'

%  A Fourier basis may also be constructed using CREATE_EASY_BASIS
%    CREATE_BASIS or MAKE_BASIS.

%  last modified 3 January 2008

%  default RANGEVAL

if nargin < 1, rangeval = [0,1];  end

%  check RANGEVAL

if length(rangeval) == 1
    if rangeval <= 0
        error('RANGEVAL a single value that is not positive.');
    end
    rangeval = [0,rangeval];
end

if rangechk(rangeval) ~= 1
    error('RANGEVAL is not a legitimate range.');
end

%  set up some default arguments

if nargin < 2, nbasis = 3;   end
if nargin < 4, dropind = []; end

%  set up period

width = rangeval(2) - rangeval(1);
if nargin < 3
    period = width;
end
if (period <= 0)
    error ('Period must be positive for a Fourier basis');
end

%  add one if NBASIS is even

if (2*floor(nbasis/2) == nbasis)
    nbasis = nbasis + 1;
end

%  check DROPIND

if length(dropind) > 0
    if length(dropind) >= nbasis
        error('Too many index values in DROPIND.');
    end
    dropind = sort(dropind);
    if length(dropind) > 1
        if any(diff(dropind)) == 0
            error('Multiple index values in DROPIND.');
        end
    end
    for i=1:length(dropind);
        if dropind(i) < 1 || dropind(i) > nbasis
            error('An index value is out of range.');
        end
    end
end

type        = 'fourier';
params      = period;
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);


