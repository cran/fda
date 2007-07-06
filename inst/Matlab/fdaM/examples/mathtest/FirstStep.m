function P0 = FirstStep(dichtest, thetaq, charwrd)
%  initialize EM algorithm by finding an initial set of probabilities
%  Arguments:
%    dichtest ... a nex by nit matrix of binary item scores
%    thetaq   ... set of quadrature points
%    charwrd  ... if 1, data are in character mode, otherwise numeric
%  Return:
%    P0 ... a nq by nit matrix of proportions 

if nargin < 3, charwrd = 1; end

[nex,nit] = size(dichtest);  % compute no. examinees and no. items
nq = length(thetaq);         %  number of trait values

% we want to make a histogram, with each bin centered
%   on a trait value.  First construct nq+1 boundaries for bars

bounds = zeros(nq+1,1);
bounds(   1) = -1e10;
bounds(nq+1) =  1e10;
for q = 2:nq
  bounds(q) = (thetaq(q-1)+thetaq(q))/2;
end

%  for each examinee, compute the index of the bar or bin
%    containing his quantile value

qscore = norminv((1:nex)./(nex+1), 0, 1);

%  get indices of bins corresponding to quantiles

binindex = zeros(nex,1);
for i=1:nq
  index = (qscore <= bounds(i+1) & qscore > bounds(i));
  binindex(index) = i;
end

%  Get preliminary estimates of IRF's at thetaq values
%  note:  the data play a role here only in terms of the
%  sorting index array, sortindex.  This array sorts the rows of
%  the dichotomous response matrix according to the rank of
%  the number right scores.

%  compute scores on the test

if charwrd
    score = zeros(nex,1);
    for i=1:nex
        temp = double(dichtest(i,:))-48;
        score(i) = sum(temp); 
    end
else
    score = sum(dichtest')'; 
end

%  sscore are the sorted scores, sortindex the indices that
%     sort vector score.  A random normal deviate, mean 0,
%     std. dev. .01 is added to each score before sorting
%     to sort tied values in random order.

[sscore,sortindex] = sort(score+0.01.*randn(nex,1));

%  for each item, compute proportion of examinees in each bin
%  that pass the test using function bin

biny = zeros(nq,nit);   
for j = 1:nit
    if charwrd
        temp = double(dichtest(sortindex,j))-48; 
    else
        temp = dichtest(sortindex,j);
    end
   biny(:,j) = bin(nq, temp, binindex);
end

%  Function bound replaces probability values in biny by 
%  1/(2*NEX) if lower, or by 1 - 1/(2*NEX) if higher.

P0 = bound(biny, nex);

