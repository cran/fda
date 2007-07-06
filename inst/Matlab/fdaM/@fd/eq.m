function fdequal = eq(fd1,fd2)
%  EQ tests for equality of two functional data objects

%  Last modfied 20 July 2006

fdequal = 1;

basis1 = getbasis(fd1);
basis2 = getbasis(fd2);

if ~(basis1 == basis2)
    fdequal = 0;
    return;
end

coef1 = getcoef(fd1);
coef2 = getcoef(fd2);

if ~isequal(coef1,coef2)
    fdequal = 0;
    return;
end

fdnames1 = getnames(fd1);
fdnames2 = getnames(fd2);

if ~strcmp(fdnames1{1},fdnames2{1}) || ...
   ~strcmp(fdnames1{2},fdnames2{2}) || ...
   ~strcmp(fdnames1{3},fdnames2{3}) 
    fdequal = 0;
    return;
end
