function basisequal = eq(basis1, basis2)

% EQ assesses whether two bases are equivalent.

%  Last modified 20 July 2006

type1   = getbasistype(basis1);
range1  = getbasisrange(basis1);
nbasis1 = getnbasis(basis1);
pars1   = getbasispar(basis1);
drop1   = getdropind(basis1);

type2   = getbasistype(basis2);
range2  = getbasisrange(basis2);
nbasis2 = getnbasis(basis2);
pars2   = getbasispar(basis2);
drop2   = getdropind(basis2);

basisequal = 1;

%  check types

if ~strcmp(type1,type2)
    basisequal = 0;
    return;
end

%  check ranges

if range1(1) ~= range2(1) || range1(2) ~= range2(2)
    basisequal = 0;
    return;
end

%  check numbers of basis functions

if nbasis1 ~= nbasis2
    basisequal = 0;
    return;
end

%  check parameter vectors

if ~isequal(pars1, pars2)
    basisequal = 0;
    return;
end

%  check indices of basis function to drop

if ~isequal(drop1, drop2)
    basisequal = 0;
    return;
end
