function [N, CN, CP, L, CL] = Estep(dichtest, P, wgtq, charwrd)
%  The E step of the EM algorithm
%  dichtest ... N by n matrix of binary item scores
%  P        ... Q by n matrix of probabilities
%  wgtq     ... Q by 1 vector of quadrature weights
%  charwrd  ... if 1 data are in character mode, otherwise numeric

if nargin < 4, charwrd = 1;  end

% get number of examinees and number of items
[nex, nit] = size(dichtest); 
Q = length(wgtq);  %  get number of quadrature weights
%  Compute:
%  CL:  N by Q matrix of conditional likelihoods, 
%  L:   N marginal likelihoods,
%  CP:  N by Q matrix of conditional probabilities, and
%  CN:  n by Q matrix of conditional pseudo-frequencies
CL = zeros(nex,Q);
CP = CL;
CN = zeros(nit,Q);
L  = zeros(nex,1);
logP   = log(P);
log1mP = log(1 - P);
for i=1:nex
    if charwrd
        temp = (double(dichtest(i,:))-48)';
    else
        temp = dichtest(i,:)';
    end
    %notmiss = (temp ~= 2);
    %temp = temp(notmiss);
    %CL(i,:) = exp(logP(:,notmiss)*temp + log1mP(:,notmiss)*(1-temp))';
    CL(i,:) = exp(logP*temp + log1mP*(1-temp))';
    %plot(thetaq,CL(i,:)');   
    L(i) = CL(i,:)*wgtq;
    CP(i,:) = wgtq'.*CL(i,:)./L(i);
    %CN(notmiss,:) = CN(notmiss,:) + temp*CP(i,:);
    CN = CN + temp*CP(i,:);
    %pause
end
%  Compute Q marginal pseudo-frequencies
N = sum(CP);

