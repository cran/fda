%  run some tests on basis functions

nbasis   = 5;

%  ------------------  Bspline basis  -----------------------

rangeval = [0,1];

clear basisobj
basisobj = create_bspline_basis()

clear basisobj
basisobj = create_bspline_basis(rangeval)

clear basisobj
basisobj = create_bspline_basis(rangeval, nbasis);
basisobj

clear basisobj
norder = 3;
basisobj = create_bspline_basis(rangeval, nbasis, norder);
basisobj

clear basisobj
norder = 4;
basisobj = create_bspline_basis(rangeval, norder+5, norder, ...
                       [0, 0.2, 0.45, 0.5, 0.55, 0.8, 1]);
basisobj

clear basisobj
basisobj = create_bspline_basis(rangeval, nbasis);

plot(basisobj)

full(eval_penalty(basisobj, 2))

clear basisobj
basisobj = create_bspline_basis(rangeval, nbasis, 4, [0,0.5,1]);
basisobj

clear basisobj
basisobj = create_bspline_basis(rangeval, nbasis, 4, [0,0.5,1], ...
                               [1,nbasis]);
basisobj

full(eval_penalty(basisobj, 2))

%  -------------------------  constant basis  --------------------

clear basisobj
basisobj = create_constant_basis();
basisobj

clear basisobj
basisobj = create_constant_basis(rangeval);
basisobj

plot(basisobj)

%  -----------------------  exponential basis  ------------------

clear basisobj
basisobj = create_exponential_basis();
basisobj

clear basisobj
basisobj = create_exponential_basis(rangeval);
basisobj

clear basisobj
basisobj = create_exponential_basis(rangeval, nbasis);
basisobj

clear basisobj
basisobj = create_exponential_basis(rangeval, nbasis, (-2:2));
basisobj

plot(basisobj)

eval_penalty(basisobj, 2)

clear basisobj
basisobj = create_exponential_basis(rangeval, nbasis, (-2:2), 3);

eval_penalty(basisobj, 2)

%  ------------------------  fourier basis  ---------------------

rangeval = [-pi,pi];

clear basisobj
basisobj = create_fourier_basis();
basisobj

clear basisobj
basisobj = create_fourier_basis(rangeval);
basisobj

clear basisobj
basisobj = create_fourier_basis(rangeval, nbasis);
basisobj

clear basisobj
basisobj = create_fourier_basis(rangeval, nbasis, pi);
basisobj

plot(basisobj)

eval_penalty(basisobj, 2)

clear basisobj
basisobj = create_fourier_basis(rangeval, nbasis, pi, 1);

eval_penalty(basisobj, 2)

%  --------------------  monomial basis  ------------------------

rangeval = [-1,1];

clear basisobj
basisobj = create_monomial_basis();
basisobj

clear basisobj
basisobj = create_monomial_basis(rangeval);
basisobj

clear basisobj
basisobj = create_monomial_basis(rangeval, nbasis);
basisobj

clear basisobj
basisobj = create_monomial_basis(rangeval, nbasis, [0,1,3,5,7]);
basisobj

plot(basisobj)

eval_penalty(basisobj, 2)

clear basisobj
basisobj = create_monomial_basis(rangeval, nbasis, ...
                                 (0:nbasis-1), 1:2);

eval_penalty(basisobj, 2)

%  -----------------------  polygonal basis  --------------------

rangeval = [0,1];

clear basisobj
basisobj = create_polygonal_basis();
basisobj

clear basisobj
basisobj = create_polygonal_basis(rangeval);
basisobj

plot(basisobj)

clear basisobj
basisobj = create_polygonal_basis(rangeval, 0.5);
basisobj

plot(basisobj)

full(eval_penalty(basisobj, 0))
full(eval_penalty(basisobj, 1))

clear basisobj
basisobj = create_polygonal_basis(rangeval, [0.25,0.5,0.75]);
basisobj

full(eval_penalty(basisobj, 0))
full(eval_penalty(basisobj, 1))

clear basisobj
basisobj = create_polygonal_basis(rangeval, [0.25,0.5,0.75], ...
                                  [1,5]);
basisobj

full(eval_penalty(basisobj, 0))
full(eval_penalty(basisobj, 1))

%  -----------------------  power basis  ------------------------

rangeval = [0.1,1];

clear basisobj
basisobj = create_power_basis();
basisobj

clear basisobj
basisobj = create_power_basis(rangeval);
basisobj

clear basisobj
basisobj = create_power_basis(rangeval, nbasis, [0, 0.5, 1, 1.5, 2]);
basisobj

plot(basisobj)

eval_penalty(basisobj, 2)

inprod(basisobj, basisobj, 2, 2)




