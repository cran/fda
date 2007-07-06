function fdParobj = putLfd(fdParobj, Lfdobj)
%  Replace the Lfd parameter and reset penmat parameter to []

%  Last modified 17 May 2004

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object');
end

fdParobj.Lfd    = Lfdobj;
fdParobj.penmat = [];