function Lfdobj = Lfd(nderiv, bwtcell)
%  LFD creates a linear differential operator object of the form
%
%  Lx(t) = w_0(t) x(t) + ... + w_{m-1}(t) D^{m-1}x(t) +
%          \exp[w_m(t) D^m x(t)  where nderiv = nderiv.
%
%  Function x(t) is operated on by this operator L, and the operator
%  computes a linear combination of the function and its first nderiv
%  derivatives.  The function x(t) must be scalar.
%
%  The linear combination of derivatives is defined by the weight
%  or coefficient functions w_j(t), and these are assumed to vary
%  over t, although of course they may also be constant as a
%  special case.
%
%  The weight coefficient for D^m is special in that it must
%  be positive to properly identify the operator.  This is why
%  it is exponentiated.  In most situations, it will be 0,
%  implying a weight of one, and this is the default.
%
%  The inner products of the linear differential operator L
%  applied to basis functions is evaluated in the functions
%  called in function EVAL_PENALTY().
%
%  Some important functions also have the capability of allowing
%  the argument that is an LFD object be an integer. They convert
%  the integer internally to an LFD object by INT2LFD().  These are:
%     EVAL_FD()
%     EVAL_MON()
%     EVAL_POS()
%     EVAL_BASIS()
%     EVAL_PENALTY()
%
%  Arguments:
%
%  NDERIV ... the order of the operator
%          the highest order of derivative.
%  BWTCELL ... A cell vector object with either NDERIV or
%              NDERIV+1 cells.
%          If there are NDERIV cells, then the coefficient of
%          D^NDERIV is set to 1; otherwise, cell NDERIV+1
%          contains a function that is exponentiated to define
%          the actual coefficient.
%
%  Simple cases:
%
%  All this generality may not be needed, and, for example,
%  often the linear differential operator will be
%  simply L = D^m, defining Lx(t) = D^mx(t).  Or the weights and
%  forcing functions may all have the same bases, in which case
%  it is simpler to use a functional data objects to define them.
%  These situations cannot be accommodated within Lfd(), but
%  there is function int2Lfd(m) that converts a nonnegative
%  integer nderiv into an Lfd object equivalent to D^m.
%  There is also fd2cell(fdobj) and that converts a functional
%  data object into cell object, which can then be used as
%  an argument of Lfd().
%
%  Returns:
%
%  LFDOBJ ... a functional data object

%  last modified 3 January 2008

%  check nderiv

if ~isnumeric(nderiv)
    error('Order of operator is not numeric.');
end
if nderiv ~= round(nderiv)
    error('Order of operator is not an integer.');
end
if nderiv < 0
    error('Order of operator is negative.');
end

%  check that BWTCELL is a cell object

if ~iscell(bwtcell)
    error('BWTCELL not a cell object.');
end

if isempty(bwtcell)
    if nderiv > 0
        error(['Positive derivative order accompanied by ', ...
               'empty weight cell.']);
    end
else
    bwtsize = size(bwtcell);
    bfdPar  = bwtcell{1};
    bfd     = getfd(bfdPar);
    brange  = getbasisrange(getbasis(bfd));

    %  BWTCELL two-dimensional.
    %  Only possibilities are (1) N > 1, nderiv == 1, and
    %                         (2) N = 1, nderiv >= 1;
    if length(bwtsize) == 2
        if bwtsize(1) > 1 &&  bwtsize(2) > 1
            error('BWTCELL is not a vector.');
        else
            if nderiv > 0
                if bwtsize(1) ~= nderiv && bwtsize(2) ~= nderiv
                    error(['Dimension of BWTCELL not ', ...
                           'compatible with NDERIV.']);
                end
            end
        end
    end

    %  BWTCELL has more than two dimensions.

    if length(bwtsize) > 2
        error('BWTCELL has more than two dimensions.');
    end

    %  Check the ranges for compatibility

    for j=2:nderiv
        brangej  = getbasisrange(getbasis(getfd(bwtcell{j})));
        if any(brangej ~= brange)
            error('Incompatible ranges in weight functions.');
        end
    end


end

%  set up the Lfd object.

Lfdobj.nderiv  = nderiv;
Lfdobj.bwtcell = bwtcell;

Lfdobj = class(Lfdobj, 'Lfd');

