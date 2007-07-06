function S = smoothmat(thetaq, N)
n = length(thetaq);
h = 1.1*N^(-.2);
phiq = zeros(1,n);
kern = zeros(n,n);
thetctr = (thetaq(2:n) + thetaq(1:(n-1)))./2;
for i=1:n
    if i == 1
        phiq(i) = normcdf(thetctr(1));
    end
    if i == n
        phiq(i) = 1 - normcdf(thetctr(n-1));
    end
    if i > 1 & i < n
        phiq(i) = normcdf(thetctr(i)) - normcdf(thetctr(i-1));
    end
    for j=1:n
        kern(i,j) = normpdf((thetaq(j) - thetaq(i))/h);
    end
end
S = zeros(n,n);
for i=1:n
    S(i,:) = phiq.*kern(i,:)./sum(phiq.*kern(i,:),2);
end
