function display(bifd)
%  DISPLAY  Display a functional data object.

%  Last modified 11 May 2004

  if strcmp(class(bifd), 'bifd')
    fprintf('Dimensions of data:\n');
    fdnames = getnames(bifd);
    fprintf(['   ',fdnames{1},'\n']);
    fprintf(['   ',fdnames{2},'\n']);
    fprintf(['   ',fdnames{3},'\n']);
  else
    error('Argument not a functional data object');
  end

  sbasisobj = bifd.sbasisobj;
  fprintf('\nBasis for first argument:\n');
  display(sbasisobj);

  tbasisobj = bifd.tbasisobj;
  fprintf('\nBasis for second argument:\n');
  display(tbasisobj);
