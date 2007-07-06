function fdnames = getnames(bifd)
%  GETNAMES Extracts the fdnames from a functional data object FD

if isa_bifd(bifd)
    fdnames = bifd.bifdnames;
else
    error('Argument is not a bivariate functional data object.');
end

