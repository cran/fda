%  This is code for testing function fSection

%  set up some data for a sine with period 1 over [0,5]

x = linspace(0,5,501)';
y = sin(2*pi*x);

%  set up a basis for these data, knot at each point,
%  order 6 so as to have a smooth 2nd derivative

basisin = create_bspline_basis([0,5],255,6);

%  smooth the data

fdin = smooth_basis(x, y, basisin);

%  set up four cutpoints to make 5 sections

cutpoints = 1:4;

%  set up order 6 bspline basis with knots at 0, 0.05, ...

basisobj = create_bspline_basis([0,1],55,6);

%  section the sine curve

fdout = fSection(fdin, cutpoints, basisobj);

% plot the sections

plot(fdout)

% plot the first derivative of the sections

plot(fdout,1)

% plot the second derivative of the sections

plot(fdout, 2)

%  set up data for the true curve for each section

argvals = linspace(0,1,101)';
sinvec  = sin(2*pi*argvals);

%  plot the fit of the first section to the true curve

plotfit_fd(sinvec, argvals, fdout(1))

%  set up data for the true second derivative

D2sinvec = -(2*pi)^2.*sinvec;

%  evaluate the second derivative of the first section

D2sinhat = eval_fd(argvals, fdout(1), 2);

%  plot the fit and true curve

plot(argvals, D2sinvec, 'b-', argvals, D2sinhat, 'r:')

%  Display the maximum error

maxerr = max(abs(D2sinhat-D2sinvec));
disp(['maximum error = ',num2str(maxerr)])

%  Display the maximum relative error

maxrelerr = maxerr/max(abs(D2sinvec));
disp(['maximum relative error = ',num2str(maxrelerr)])


