function display(basis)
%  DISPLAY  Display a functional data basis object.

%  Last modified 22 December 2007

if ~strcmp(class(basis), 'basis')
    error('Argument not a functional data object');
end

fprintf('\nBasis:\n');

%  display type

fprintf(['  Type: ', basis.type,'\n']);

%  display range

fprintf(['  Range: ', num2str(basis.rangeval(1)), ...
         ' to ',      num2str(basis.rangeval(2)),'\n']);

%  return if a constant basis

if strcmp(basis.type, 'const')
    return;
end

%  display number of basis functions

fprintf(['  Number of basis functions: ', ...
         num2str(basis.nbasis),     '\n']);

%  display parameters according to type of basis

if  strcmp(basis.type, 'fourier')
    fprintf(['  Period: ',num2str(basis.params),'\n']);
end
if strcmp(basis.type, 'bspline')
    norder = basis.nbasis - length(basis.params);
    fprintf(['  Order of spline: ', ...
            num2str(norder),     '\n']);
    if length(basis.params) > 0
        fprintf('  Interior knots\n');
        disp(basis.params);
    else
        fprintf('  There are no interior knots.\n');
    end
end
if strcmp(basis.type, 'polyg')
    fprintf('  Argument values\n');
    disp(basis.params);
end
if strcmp(basis.type, 'expon')
    fprintf('  Rate coefficients\n');
    disp(basis.params);
end
if strcmp(basis.type, 'power')
    fprintf('  Exponents\n');
    disp(basis.params);
end

%  display indices of basis functions to be dropped

if length(basis.dropind) > 0
    fprintf('  Indices of basis functions to be dropped\n');
    disp(basis.dropind);
end
