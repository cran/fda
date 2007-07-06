function isafdPar = isa_fdPar(fdParstr)
%  ISA_FDPAR  checks an object to see if it is of the FD class 

%  last modified 9 September 2003

  isafdPar = 1;
  if ~(strcmp(class(fdParstr), 'fdPar') | ...
       strcmp(class(fdParstr), 'fdpar'))
    isafdPar = 0;
    return;
  end

