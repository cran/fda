function Xmat = symsolve(Asym, Bmat)
%SYMSOLVE solves ASYM Xmat = BMAT for Xmat where ASYM is symmetric
%  Last modified 24 April 2003

if max(max(abs(Asym-Asym')))/max(max(abs(Asym))) > 1e-10
    error('First argument is not symmetric.');
else
    Asym = (Asym + Asym')./2;
end
%  Choleski decomposition Asym = Rmat'*Rmat;
[Rmat,p] = chol(Asym);
if p > 0
    error('Argument ASYM is singular.');
end
temp = Rmat'\Bmat;
Xmat = Rmat\temp;
