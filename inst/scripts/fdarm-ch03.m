%
%
% Ramsay, Hooker & Graves (2009)
% Functional Data Analysis with R and Matlab (Springer)
%

%  Remarks and disclaimers

%  These R commands are either those in this book, or designed to 
%  otherwise illustrate how R can be used in the analysis of functional
%  data.  
%  We do not claim to reproduce the results in the book exactly by these 
%  commands for various reasons, including:
%    -- the analyses used to produce the book may not have been
%       entirely correct, possibly due to coding and accuracy issues
%       in the functions themselves 
%    -- we may have changed our minds about how these analyses should be 
%       done since, and we want to suggest better ways
%    -- the R language changes with each release of the base system, and
%       certainly the functional data analysis functions change as well
%    -- we might choose to offer new analyses from time to time by 
%       augmenting those in the book
%    -- many illustrations in the book were produced using Matlab, which
%       inevitably can imply slightly different results and graphical
%       displays
%    -- we may have changed our minds about variable names.  For example,
%       we now prefer "yearRng" to "dayrange" for the weather data.
%    -- three of us wrote the book, and the person preparing these scripts
%       might not be the person who wrote the text
%  Moreover, we expect to augment and modify these command scripts from time
%  to time as we get new data illustrating new things, add functionality
%  to the package, or just for fun.

%
% ch. 3.  How to specify basis systems for building functions
%

%  Set up some strings for constructing paths to folders.
%  These strings should be modified so as to provided access
%  to the specified folders on your computer.

%  Path to the folder containing the Matlab functional data analysis
%  software

fdaMPath = 'c:/Program Files/MATLAB/R2009a/fdaM';

addpath(fdaMPath)

%  Path to the folder containing the examples

examplesPath = [fdaMPath,'/examples'];

addpath(examplesPath)

%
% Section 3.1 Basis Function Systems for Constructing Functions
%

unitRng = [0,1];

const_basis   = create_constant_basis(unitRng);

nbasis = 1;
monom_basis   = create_monomial_basis(unitRng, nbasis);

nbasis = 5;
period = 1;
fourier_basis = create_fourier_basis(unitRng, nbasis, period);

nbasis = 5;
norder = 2;
breaks = linspace(0, 1, 5);
bspline_basis = create_bspline_basis(unitRng, nbasis, norder, breaks);

%
% Section 3.2 Fourier Series for Periodic Data and Functions
%

yearRng = [0,365];
nbasis  = 65;
daybasis65 = create_fourier_basis(yearRng, nbasis);

nbasis = 3;
period = 500;
daybasis_T = create_fourier_basis(yearRng, nbasis, period);

dropind = 1;
nbasis  = 65;
zerobasis1 = create_fourier_basis(yearRng, nbasis, dropind);
zerobasis2 = daybasis65(2:65);

nbasis = 3;
fourier_basis3 = create_fourier_basis(unitRng, nbasis);

help(create_fourier_basis)

%
% Section 3.3 Spline series for Non-periodic Data and Functions
%

%  section 3.3.3 Examples

%  order 4 spline, one interior knot

breaks = [0, .5, 1];
norder = 4;
nbasis = norder + length(breaks) - 2;
bspline4 = create_bspline_basis(unitRng, nbasis, norder, breaks);

plot(bspline4)

%  order 2 spline, one interiot knot

breaks = [0, .5, 1];
norder = 2;
nbasis = norder + length(breaks) - 2;
bspline2 = create_bspline_basis(unitRng, nbasis, norder, breaks);

plot(bspline2)

%  order 2 spline,  2 equal interior knots

breaks = [0, .5, .5, 1];
norder = 2;
nbasis = norder + length(breaks) - 2;
bspline2_2 = create_bspline_basis(unitRng, nbasis, norder, breaks);

plot(bspline2_2)

%  order 4 spline, 3 equal interior knots

breaks = [0, .5, .5, .5, 1];
norder = 2;
nbasis = norder + length(breaks) - 2;
bspline4_3 = create_bspline_basis(unitRng, nbasis, norder, breaks);

plot(bspline4_3)

%  section 3.3.4 B-Splines

splinebasis = create_bspline_basis([0,10], 13);

% Figure 3.1

plot(splinebasis) 

% Figure 3.2

twopiRng = [0,2*pi];
basis2 = create_bspline_basis(twopiRng, 5, 2);
basis3 = create_bspline_basis(twopiRng, 6, 3);
basis4 = create_bspline_basis(twopiRng, 7, 4);

theta     = linspace(0, 2*pi, 201);
sin_theta = sin(theta);

sin2 = smooth_basis(theta, sin_theta, basis2);
sin3 = smooth_basis(theta, sin_theta, basis3);
sin4 = smooth_basis(theta, sin_theta, basis4);

sin2_theta = eval_fd(sin2, theta);
sin3_theta = eval_fd(sin3, theta);
sin4_theta = eval_fd(sin4, theta);

sinRng = [min(sin2_theta), max(sin2_theta)];
pi3    = (1:3)*pi/2;

subplot(3,2,1)
plot(theta, sin2_theta, 'b-', theta, sin_theta, 'b--')
hold on
for i=1:3
    plot([pi3(i), pi3(i)], sinRng, 'b:') 
end
hold off
xlabel('') 
ylabel('Order = 2')
title('sine(t)')
axis([twopiRng,sinRng])

Dsin2_theta = eval_fd(sin2, theta, 1);
subplot(3,2,2)
plot(theta, Dsin2_theta, 'b-', theta, cos(theta), 'b--')
hold on
for i=1:3
    plot([pi3(i), pi3(i)], sinRng, 'b:') 
end
hold off
xlabel('') 
ylabel('Order = 2')
title('D sine(t)')
axis([twopiRng,sinRng])

subplot(3,2,3)
plot(theta, sin3_theta, 'b-', theta, sin_theta, 'b--')
hold on
for i=1:3
    plot([pi3(i), pi3(i)], sinRng, 'b:') 
end
hold off
xlabel('') 
ylabel('Order = 3')
axis([twopiRng,sinRng])

Dsin3_theta = eval_fd(sin3, theta, 1);
subplot(3,2,4)
plot(theta, Dsin3_theta, 'b-', theta, cos(theta), 'b--')
hold on
for i=1:3
    plot([pi3(i), pi3(i)], sinRng, 'b:') 
end
hold off
xlabel('') 
ylabel('Order = 3')
axis([twopiRng,sinRng])

subplot(3,2,5)
plot(theta, sin4_theta, 'b-', theta, sin_theta, 'b--')
hold on
for i=1:3
    plot([pi3(i), pi3(i)], sinRng, 'b:') 
end
hold off
xlabel('') 
ylabel('Order = 4')
axis([twopiRng,sinRng])

Dsin4_theta = eval_fd(sin4, theta, 1);
subplot(3,2,6)
plot(theta, Dsin4_theta, 'b-', theta, cos(theta), 'b--')
hold on
for i=1:3
    plot([pi3(i), pi3(i)], sinRng, 'b:') 
end
hold off
xlabel('') 
ylabel('Order = 4')
axis([twopiRng,sinRng])

%  order 6 spline with equally spaced knots

splinebasis6 = create_bspline_basis([0,10], 15, 6);

%  plot a central basis function

subplot(1,1,1)
plot(splinebasis6(7))

%  plot two left boundary basis functions

plot(splinebasis6(1:2))

%  order 4 spline with knots not equally spaced

breaks=[0, .7, 1];
norder = 4;
nbasis = norder + length(breaks) - 2;
spline_unequal_knots = ...
    create_bspline_basis(unitRng, nbasis, norder, breaks);

%  plot sum of basis function values (all equal to 1)

t10 = 0:0.1:10;
spline6_10 = eval_basis(t10, splinebasis6);
plot(t10, sum(spline6_10,2))

%
% Section 3.4 Constant, Monomial and Other Bases
%

conbasis = create_constant_basis(unitRng);
conb     = create_constant_basis();
disp(eq(conbasis, conb))

nbasis = 4;
monbasis = create_monomial_basis(unitRng, nbasis);

%
% Section 3.5 Methods for Functional Basis Objects
%

disp(monbasis)

disp(monbasis == monbasis)

isa_basis(monbasis)

disp(getbasispar(monbasis))

t_1        = 0:0.1:1;
monbmat_1  = eval_basis(t_1, monbasis);
Dmonbmat_1 = eval_basis(t_1, monbasis, 1);

%
% Section 3.6 The Structure of the basisfd or basis Class
%

help(basis)

%
% Section 3.7 Some Things to Try
%
% (exercises for the reader)
