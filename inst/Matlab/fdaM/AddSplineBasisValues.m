function newbasisobj = AddSplineBasisValues(basisobj, nquad, tobs, nderiv)

% of quadrature points spanning an interval between knots.

if nargin < 4,  nderiv = 0;   end
if nargin < 3,  tobs   = [];  end
if nargin < 2,  nquad  = 0;   end

type = getbasistype(basisobj);
if ~strcmp(type, 'bspline')
    error('BASISOBJ is not a bspline basis object.');
end

newbasisobj = basisobj;

rangeval = getbasisrange(basisobj);
params   = getbasispar(basisobj);
breaks   = [rangeval(1), unique(params), rangeval(2)];
nbreaks  = length(breaks);

%  set up quadrature points and weights using Simpson's rule

if nquad > 0
    quadpts1 = linspace(breaks(1),breaks(2),nquad)';
    quadwts1 = ones(nquad,1);
    quadwts1(2:2:nquad-1) = 4;
    quadwts1(3:2:nquad-2) = 2;
    quadwts1 = ((breaks(2)-breaks(1))/(nquad-1)).*quadwts1/3;
    quadvals = [quadpts1, quadwts1];
    
    for i=3:nbreaks
        quadptsi = linspace(breaks(i-1),breaks(i),nquad)';
        quadwtsi = ones(nquad,1);
        quadwtsi(3:2:nquad-2) = 2;
        quadwtsi(2:2:nquad-1) = 4;
        quadwtsi = ((breaks(i)-breaks(i-1))/(nquad-1)).*quadwtsi/3;
        quadvals = [quadvals;[quadptsi, quadwtsi]];
    end    
    newbasisobj = putquadvals(newbasisobj, quadvals);
end


if isempty(tobs) && nquad > 0
    %  set up values cell array for values at quadrature points only
    quadpts = quadvals(:,1);
    values = cell(1,nderiv+1);
    for i=0:nderiv
        values{i+1} = eval_basis(quadpts, basisobj, i);
    end
    newbasisobj = putvalues(newbasisobj, values);
elseif ~isempty(tobs) && nquad == 0
    basisvalues = cell(1,nderiv+2);
    % Basis function values at observation times
    basisvalues{1} = tobs;
    for i=0:nderiv
        basisvalues{i+2} = eval_basis(tobs, basisobj, i);
    end
    newbasisobj = putbasisvalues(newbasisobj, basisvalues);
else
    quadpts = quadvals(:,1);
    basisvalues = cell(2,nderiv+2);
    % Basis function values at observation times
    basisvalues{1,1} = tobs;
    for i=0:nderiv
        basisvalues{1,i+2} = eval_basis(tobs, basisobj, i);
    end
    % Basis function values at quadrature points
    basisvalues{2,1} = quadpts;
    for i=0:nderiv
        basisvalues{2,i+2} = eval_basis(quadpts, basisobj, i);
    end
    % Store the cell array in the basis object
    newbasisobj = putbasisvalues(newbasisobj, basisvalues);
end