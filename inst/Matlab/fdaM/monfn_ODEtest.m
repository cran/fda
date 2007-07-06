x = linspace(0,1,51)';
y = exp(x) - 1;

nbasis = 4;
wbasis = create_bspline_basis([0,1],nbasis);

Wfd = fd(ones(nbasis,1),wbasis);

z = monfn_ODE(x, Wfd);

