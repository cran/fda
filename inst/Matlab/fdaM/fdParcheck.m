function fdParobj = fdParcheck(fdParobj)

if ~isa_fdPar(fdParobj) 
    if isa_fd(fdParobj) || isa_basis(fdParobj)
        fdParobj = fdPar(fdParobj);
    else
        error(['FDPAROBJ is not a functional parameter object, ', ...
               'not a functional data object, and ', ...
               'not a basis object.']);
    end
end
