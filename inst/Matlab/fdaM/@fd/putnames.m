function fd = putnames(fd, fdnames)
% PUTNAMES  Assigns fdnames to a functional data object FD

%  Last modified 20 July 2006

if isa_fd(fd)
    if iscell(fdnames)
        fd.fdnames = fdnames;
    else
        error('Argument FDNAMES is not a cell object.');
    end
else
    error('Argument FD is not a functional data object.');
end

