function display(fd)
%  DISPLAY  Display a functional data object.

%  Last modified 11 May 2004

if strcmp(class(fd), 'fd')
    fprintf('Dimensions of data:\n');
    fdnames = getnames(fd);
    fprintf(['   ',fdnames{1},'\n']);
    if iscell(fdnames{2})
        fprintf(['   ',fdnames{2}{1},'\n']);
    else
        fprintf(['   ',fdnames{2},'\n']);
    end
    if iscell(fdnames{3})
        fprintf(['   ',fdnames{3}{1},'\n']);
    else
        fprintf(['   ',fdnames{3},'\n']);
    end
else
    error('Argument not a functional data object');
end

basisobj = getbasis(fd);
display(basisobj);
