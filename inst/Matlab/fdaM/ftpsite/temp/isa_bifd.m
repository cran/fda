function isabifd = isa_bifd(bifd)
%  ISA_BIFD  checks a struct object for fields 'coef' and 'basisstr'

%  last modified 1 July 1998

  isabifd = 1;
  if ~strcmp(class(bifd), 'bifd')
    isabifd = 0;
    return;
  end
