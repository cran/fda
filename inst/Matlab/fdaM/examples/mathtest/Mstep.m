function newcoef = Mstep(CN, N, P, coef, phimat, penmat)
%  M-step
%  for each item in turn, minimize criterion Fj with respect to
%    values of coefficients cj
%  CN     ... nit by Q matrix of conditional pseudo frequencies
%  N      ... vector of Q marginal pseudo frequencies
%  P      ... Q by nit matrix of probabilities
%  coef   ... K by nit matrix of coefficients for basis functions
%  phimat ... Q by K matrix of basis function values
%  penmat ... K by K matrix for penalizing coefficient roughness
nit = size(P,2);
newcoef = coef;
iterlim = 10;
for j=1:nit
    cj   = coef(:,j);             % initial coefficients
    pj   = P(:,j); qj = 1-P(:,j); % initial probabilities    
    wj   = diag(pj.*qj.*N');      % initial weights
    CNj  = CN(j,:)';              % "Data" for this item
    resj = CNj - pj.*N';          % initial residuals
    %  initial function value
    Fj = -sum(CNj.*log(pj)+(N'-CNj).*log(qj)) + cj'*penmat*cj;
    %  gradient vector
    Gmat = -phimat' * resj        + 2.*penmat*cj;
    %  Hessian matrix
    Hmat =  phimat' * wj * phimat + 2.*penmat;
    %  search direction
    delta = Hmat \ Gmat;
    %  initial step size along direction
    alpha = 1; 
    Fjnew = Fj + 1;
    %  take step, and if new function value less than old, quit
    %  otherwise halve step and try again.
    %  in any case, stop when step size gets too small
    iter = 0;
    while (Fjnew > Fj & iter <= iterlim)
       iter = iter + 1;
       cjnew = cj - alpha .* delta; % update coefficients
       %  Don't let coefficients get too large in either direction
       cjnew(cjnew > 20) = 20; cjnew(cjnew < -20) = -20;
       % update probabilities
       pjnew = 1./(1+exp(-phimat*cjnew)); qjnew = 1 - pjnew;
       % compute new function value
       Fjnew = -sum(CNj.*log(pjnew)+(N'-CNj).*log(qjnew)) + ...
                 cjnew'*penmat*cjnew;
       %fprintf('%g ', [j, iter, Fjnew, Fj]); fprintf('\n');
       alpha = alpha/2; % halve step size
    end
    newcoef(:,j) = cjnew;   % replace old coefficients by new ones
end
