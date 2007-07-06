function display(fd)
%  DISPLAY  Display a functional data object.

%  Last modified 11 May 2004

if strcmp(class(fd), 'fd')
    fprintf('Dimensions of data:\n');
    fdnames = getnames(fd);
    fprintf(['   ',fdnames{1},'\n']);
    fprintf(['   ',fdnames{2},'\n']);
    fprintf(['   ',fdnames{3},'\n']);
else
    error('Argument not a functional data object');
end

basisobj = getbasis(fd);
display(basisobj);
