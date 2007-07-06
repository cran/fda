function display(bifdPar)
% displays bivariate functional parameter object

%  last modified 14 October 2003

fprintf('\nFD:\n');
display(bifdPar.bifd);

fprintf('\nestimate = %d\n', bifdPar.estimate);
display(bifdPar.estimate);

fprintf('\nlambdas = %.6g\n');
display(bifdPar.lambdas);
fprintf('\nlambdat = %.6g\n');
display(bifdPar.lambdat);

fprintf('\nLfds:\n');
display(bifdPar.Lfds);
fprintf('\nLfdt:\n');
display(bifdPar.Lfdt);

