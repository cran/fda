function isafd = isa_fd(fdstr)
%  ISA_FD checks a struct object for having the FD class

%  last modified 9 Septebmer 2003

  isafd = 1;
  if ~strcmp(class(fdstr), 'fd')
    isafd = 0;
    return;
  end

