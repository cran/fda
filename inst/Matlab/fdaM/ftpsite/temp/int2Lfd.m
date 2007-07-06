function Lfdobj = int2Lfd(m)
%INT2LFD converts a nonnegative integer to a linear differential
%  operator object that is equivalent to D^m.  The range of the
%  functional data object in any cell is set to [0,1], and is
%  not actually used when a linear differential operator object
%  of this nature is applied.  
%  In the event that m is already a linear differential operator
%  object, it returns the object immediately.  Thus, INT2LFD can
%  be used to screen whether an object is an integer or not.

%  Last modified 10 May 2004

%  check M

if isa_Lfd(m)
    Lfdobj = m;
    return;
end

if ~isnumeric(m)
    error(['Argument not numeric ', ...
           'and not a linear differential operator.']);
end

if length(m) ~= 1
    error('Argument is not a scalar.');
end
if round(m) ~= m
    error('Argument is not an integer.');
end
if m < 0
    error('Argument is negative.');
end

%  all the checks passed, set up a functional data object
%  The range will not be used in this case, and can be set
%  to [0, 1]

%  set up the cell object for the homogeneous part

if m == 0
    %  if derivative is zero, WFDCELL is empty
    wfdcell = {};
else
    wfd0 = fd(0,create_constant_basis([0,1]));
    wfdcell{1} = fdPar(wfd0);
    for j=2:m
        wfdcell{j} = fdPar(wfd0);
    end
end

%  define the Lfd object

Lfdobj = Lfd(m, wfdcell);
