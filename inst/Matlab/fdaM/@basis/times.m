function prodbasisobj = times(basisobj1, basisobj2)
% TIMES for two basis objects sets up a basis suitable for 
%  expanding the pointwise product of two functional data
%  objects with these respective bases.  
% In the absence of a true product basis system in this code,
%  the rules followed are inevitably a compromise:
%  (1) if both bases are B-splines, the norder is the sum of the
%      two orders - 1, and the breaks are the union of the
%      two knot sequences, each knot multiplicity being the maximum
%      of the multiplicities of the value in the two break sequences.
%      That is, no knot in the product knot sequence will have a
%      multiplicity greater than the multiplicities of this value
%      in the two knot sequences.  
%      The rationale this rule is that order of differentiability
%      of the product at each value will be controlled  by
%      whichever knot sequence has the greater multiplicity.  
%      In the case where one of the splines is order 1, or a step
%      function, the problem is dealt with by replacing the
%      original knot values by multiple values at that location
%      to give a discontinuous derivative.
%  (2) if both bases are Fourier bases, AND the periods are the 
%      the same, the product is a Fourier basis with number of
%      basis functions the sum of the two numbers of basis fns.
%  (3) if only one of the bases is B-spline, the product basis
%      is B-spline with the same knot sequence and order two
%      higher.
%  (4) in all other cases, the product is a B-spline basis with
%      number of basis functions equal to the sum of the two
%      numbers of bases and equally spaced knots.  

%  Of course the ranges must also match.

%  Last modified 22 March 2007

%  check the ranges

range1 = getbasisrange(basisobj1);
range2 = getbasisrange(basisobj2);
if range1(1) ~= range2(1) || range1(2) ~= range2(2)
    error('Ranges are not equal.');
end

%  get the types

type1 = getbasistype(basisobj1);
type2 = getbasistype(basisobj2);

%  deal with constant bases

if strcmp(type1, 'const') && strcmp(type2, 'const')
    prodbasisobj = create_constant_basis(range1);
    return;
end

if strcmp(type1, 'const')
    prodbasisobj = basisobj2;
    return;
end

if strcmp(type2, 'const')
    prodbasisobj = basisobj1;
    return;
end

%  get the numbers of basis functions including dropped indices

nbasis1 = getnbasis(basisobj1) + length(getdropind(basisobj1));
nbasis2 = getnbasis(basisobj2) + length(getdropind(basisobj2));

%  work through the cases

if strcmp(type1, 'bspline') && strcmp(type2, 'bspline')
    %  both bases B-splines
    %  get orders
    interiorknots1 = getbasispar(basisobj1);
    interiorknots2 = getbasispar(basisobj2);
    uniqueknots = union(interiorknots1,interiorknots2);
    nunique = length(uniqueknots);
    multunique = zeros(nunique,1);
    for i=1:nunique
        mult1 = length(find(interiorknots1==uniqueknots(i)));
        mult2 = length(find(interiorknots2==uniqueknots(i)));
        multunique(i) = max(mult1,mult2);
    end
    allknots = zeros(sum(multunique),1);
    m2 = 0;
    for i=1:nunique
        m1 = m2 + 1;
        m2 = m2 + multunique(i);
        allknots(m1:m2) = uniqueknots(i);
    end    
    norder1 = nbasis1 - length(interiorknots1);
    norder2 = nbasis2 - length(interiorknots2);
    norder  = norder1 + norder2 - 1;
    breaks  = [range1(1), allknots', range1(2)];
    nbasis  = length(breaks) + norder - 2;
    prodbasisobj = ...
        create_bspline_basis(range1, nbasis, norder, breaks);
    return;
end

if strcmp(type1, 'fourier') && strcmp(type2, 'fourier')
    %  both bases Fourier
    %  check whether periods match
    %  if they do not, default to the basis below.
    period1 = getbasispar(basisobj1);
    period2 = getbasispar(basisobj2);
    nbasis  = nbasis1 + nbasis2;
    if period1 == period2
        prodbasisobj = ...
            create_fourier_basis(range1, nbasis, period1);
        return;
    end
end

%  default case when all else fails: the product basis is B-spline 
%  When neither basis is a B-spline basis, the order
%  is the sum of numbers of bases, but no more than 8.
%  When one of the bases if B-spline and the other isn't,
%  the order is the smaller of 8 or the order of the spline
%  plus 2.

if strcmp(type1, 'bspline') || strcmp(type2, 'bspline')
    norder = 8;
    if strcmp(type1, 'bspline') 
        interiorknots1 = getbasispar(basisobj1);
        norder1 = nbasis1 - length(interiorknots1);
        norder = min(norder1+2, norder);
    end
    if strcmp(type2, 'bspline') 
        interiorknots2 = getbasispar(basisobj2);
        norder2 = nbasis2 - length(interiorknots2);
        norder = min(norder2+2, norder);
    end
else
    %  neither basis is B-spline
    norder = min(8, nbasis1+nbasis2);
end
%  set up the default B-spline product basis
nbasis = max(nbasis1+nbasis2, norder+1);
prodbasisobj = create_bspline_basis(range1, nbasis, norder);