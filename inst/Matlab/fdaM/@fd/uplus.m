function plusfd = uplus(fdobj)
% Unary plus of functional data object.

%  last modified 20 July 2006

if ~isa_fd(fdobj)
    error('Argument is not a functional data object.');
else
    plusfd  = fdobj;
end


