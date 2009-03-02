function newbasisobj = putquadvals(basisobj, quadvals)
%  PUQUADVALS   Enters the quadrature points and weights
%     from basis object BASISOBJ into slot basisobj.quadvals
%
%  Arguments:
%  BASISOBJ ... A basis object
%  QUADVALS ... Either a NQUAD by 2 matrix, the first column
%               containing quadrature points and the second
%               quadrature weights, 
%               or an integer greater than 1, in qudrature
%               points and weights are set to Simpson's rule.
%               Do not use this default option, however, if
%               there are repeated knot values.  

%  last modified 18 June 2007

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

%  check quadvals

if nargin < 2 || isempty(quadvals)
    %  set up default quadrature values using 20 times
    %  number of intervals
    nquad = (length(basisobj.params)+1)*20 + 1;
    if 2*floor(nquad/2) == nquad
        nquad = nquad + 1;
    end
    quadvals = genquadvals(basisobj, nquad);
else
    [nrow,ncol] = size(quadvals);
    if nrow == 2 && ncol > 2
        quadvals = quadvals';
    end
    nquad = size(quadvals,1);
    if nrow == 1 && ncol == 1
        if floor(quadvals) == quadvals && quadvals > 4
            nquad = quadvals;
            if 2*floor(nquad/2) == nquad
                nquad = nquad + 1;
            end
            quadvals = genquadvals(basisobj, nquad);
        else
            error(['A single number supplied for QUADVALS that is ', ...
                   'not a correct number of quadrature values.']);
        end
    end
    if nquad < 2
        error('Less than two quadrature points are supplied.');
    end
    if ncol ~= 2
        error('QUADVALS does not have two columns.');
    end
end

newbasisobj.type        = basisobj.type;
newbasisobj.rangeval    = basisobj.rangeval;
newbasisobj.nbasis      = basisobj.nbasis;
newbasisobj.params      = basisobj.params;
newbasisobj.dropind     = basisobj.dropind;
newbasisobj.quadvals    = quadvals;
newbasisobj.values      = basisobj.values;
newbasisobj.basisvalues = basisobj.basisvalues;

newbasisobj = class(newbasisobj, 'basis');

%  ----------------------------------------------------------------

function quadvals = genquadvals(basisobj, nquad)
%  generates quadrature points and values according to
%  Simpson's rule

%  check NQUAD for being even and at least 5
if nquad < 5 || 2*floor(nquad/2) == nquad
    error(['Simpsons rule quadrature points and weights cannot ', ...
          'be generated for even NQUAD or NQUAD < 5.']);
end
%  proceed to generate quadrature points and weights          
Twide = basisobj.rangeval(2)-basisobj.rangeval(1);
delta = Twide/(nquad-1);
quadvals = zeros(nquad,2);
quadvals(:,1) = linspace(0,Twide, nquad)';
quadvals(:,2) = ones(nquad,1);
quadvals(2:2:nquad-1,2) = 4;
quadvals(3:2:nquad-2,2) = 2;
quadvals(:,2) = delta.*quadvals(:,2)/3;
