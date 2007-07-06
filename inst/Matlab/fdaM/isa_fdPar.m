function isafdPar = isa_fdPar(fdParstr)
%  ISA_FDPAR  checks an object to see if it is of the FD class 

%  last modified 2 December 2006

  isafdPar = 1;
  if ~(strcmp(class(fdParstr), 'fdPar') || ...
       strcmp(class(fdParstr), 'fdpar'))
    isafdPar = 0;
    return;
  end

