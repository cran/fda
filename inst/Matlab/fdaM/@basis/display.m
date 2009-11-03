function display(basis)
%  DISPLAY  Display a functional data basis object.

%  Last modified 25 May 2010

if ~strcmp(class(basis), 'basis')
    error('Argument not a functional data object');
end

fprintf('\nBasis:\n');

%  display type

fprintf(['  Type: ', basis.type,'\n']);

%  display range

if ~isempty(basis.rangeval)
    fprintf(['  Range: ', num2str(basis.rangeval(1)), ...
        ' to ',      num2str(basis.rangeval(2)),'\n']);
end

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
if strcmp(basis.type, 'slide')
    breaks = [basis.rangeval(1), basis.params(1:(basis.nbasis-1)), ...
              basis.rangeval(2)];
    rates  = basis.params(basis.nbasis:(2*basis.nbasis-1));
    fprintf('   Breaks\n');
    disp(breaks)
    fprintf('   Rates\n');
    disp(rates)
end
if strcmp(basis.type, 'fd')
    fprintf('  Functional data object\n');
    disp(basis.params);
end
if strcmp(basis.type, 'FEM')
    fprintf('  FEM\n');
    FEMstruct = basis.params;
    disp('Points:')
    disp(FEMstruct.p);
    disp('Edges:')
    disp(FEMstruct.e);
    disp('Triangles:')
    disp(FEMstruct.t);
end

%  display indices of basis functions to be dropped

if length(basis.dropind) > 0
    fprintf('  Indices of basis functions to be dropped\n');
    disp(basis.dropind);
end
