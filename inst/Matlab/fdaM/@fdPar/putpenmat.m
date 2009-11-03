function newfdParobj = putpenmat(fdParobj, penmat)
%  Replace the PENMAT parameter

%  last modified 17 May 2003

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object');
end

newfdParobj = fdPar(fdParobj.fd,       ...
                    fdParobj.Lfd,      ...
                    fdParobj.lambda,   ...
                    fdParobj.estimate, ...
                    penmat);

