function display(fdParobj)

%  Last modified 3 January 2008

fprintf('\nFunctional parameter object:\n\n')
fprintf('\nFunctional data object:\n');
display(fdParobj.fd);
nderiv = getnderiv(fdParobj.Lfd);
if nderiv > 0
    fprintf('\nLinear differential operator object:\n\n');
    display(fdParobj.Lfd);
else
    fprintf('\nLFD      = 0');
end
fprintf('\n\nSmoothing parameter   = %.6g', fdParobj.lambda);
fprintf('\nEstimation status = %d',   fdParobj.estimate);
if ~isempty(fdParobj.penmat)
    fprintf('\n\nPenalty matrixX\n');
    disp(fdParobj.penmat)
end

