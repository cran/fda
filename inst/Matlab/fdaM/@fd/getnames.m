function fdnames = getnames(fd)
%  GETNAMES Extracts the fdnames from a functional data object FD

  if isa_fd(fd)
    fdnames = fd.fdnames;
  else
    error('Argument is not a functional data object.');
  end

