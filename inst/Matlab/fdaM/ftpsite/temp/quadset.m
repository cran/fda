function [basisobj, quadpts, quadwts] = quadset(nquad, basisobj, breaks)

%  Set up quadrature points and weights for Simpson's rule
%  and store information in BASISOBJ.  Simpson's rule is used to integrate
%  a function between successive values in vector BREAKS.
%  That is, over each interval [BREAKS(I),BREAKS(I+1)]. 
%  Simpson's rule uses NQUAD equally spaced quadrature points over this 
%  interval, starting with the the left boundary and ending with the right 
%  boundary.  The quadrature weights are the values 
%     DELTA.*[1,4,2,4,2,4,..., 2,4,1] where DELTA is the difference 
%     between successive quadrature points, that is,
%     DELTA = (BREAKS(I+1)-BREAKS(I))/(NQUAD-1).  
%  For example, if NQUAD = 7, and BREAKS = [0, .3, 1], then
%  over [0, 0.3] the quadrature points are the 7 values
%      [0, 0.05, 0.01, 0.15, 0.20, 0.25, 0.3]
%  and the quadrature weights are
%      0.05.*[1, 4, 2, 4, 2, 4, 1]
%  and defined similarly for the interval [0.3, 1], except that
%  the spacing is now 0.7/6.  

%  Arguments:
%  NQUAD    ... The number of Simpson's rule quadrature points
%               over any intervals [breaks(i),breaks(i+1)].
%               This must be at least 5 and must be odd.
%  BASISOBJ ... The basis object that will contain the 
%               quadrature points and weights
%  BREAKS   ... If this is input, these are the interval boundaries.
%               The first value must be the initial point of the
%               interval over which the basis is defined, and the
%               final value must be the end point. 
%               If this is not supplied, and BASISOBJ is of the
%               'bspline' type, then the knots are used as these
%               values.

%  Last modified 14 August 2006

%  check BASISOBJ

if ~isa_basis(basisobj)
    error('BASISOBJ is not a basis object.');
end

if nargin < 3
    type     = getbasistype(basisobj);
    if ~strcmp(type, 'bspline')
        error('BREAKS not supplied and BASISOBJ is not a spline basis.');
    end
    rangeval = getbasisrange(basisobj);
    params   = getbasispar(basisobj);
    knots    = [rangeval(1), params, rangeval(2)];
    breaks   = unique(knots);
end

nbreaks = length(breaks);

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
quadpts = quadvals(:,1);
quadwts = quadvals(:,2);

basisobj = putquadvals(basisobj, quadvals);

values = cell(2,1);
for ivalue=1:2
    basisvalues    = eval_basis(quadpts, basisobj, ivalue-1);
    values{ivalue} = basisvalues;
end

basisobj = putvalues(basisobj, values);

