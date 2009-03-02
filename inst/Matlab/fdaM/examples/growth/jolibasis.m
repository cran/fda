addpath ('c:\matlab\fdaM')
addpath ('c:\matlab\growmatlab')

%  The Jolicoeur model is:
%  h(t) = a*[1 - 1/F(t)] where
%    F(t) = [b_1*(t + e)]^c_1 + [b_2*(t + e)]^c_2 + [b_3*(t + e)]^c_3

%  the coefficients are in the order:

%          a, b_1, b_2, b_3, c_1, c_2, c_3, e

fid      = fopen('FELSJPAF.PAR','rt');
tempvec  = fscanf(fid,'%f');
ncaseJ   = length(tempvec)/21;
Felsjpaf = reshape(tempvec, [21, ncaseJ]);
clear tempvec;
parmatf = Felsjpaf(7:14,:)';

%  average coefficients for females:

Coefmean = [ ...
 164.7077, 0.3071, 0.1106, 0.0816, 0.7285, 3.6833, 16.6654, 1.4744]';

a = Coefmean(1);
b = Coefmean(2:4);
c = Coefmean(5:7);
e = Coefmean(8);

age = (0:0.1:18)';
n = length(age);

ageper = age + e;

[hgtmean, vel, acc, Fvec, fmat] = jolifn(age, Coefmean);
plot(age, hgtmean)

fmat2 = ones(n,3);

for i=1:3
    fmat2(:,i) = (b(i).*ageper).^c(i);
end

Fvec2 = sum(fmat2,2);

fmat2 = [fmat2, -ones(n,1)];

denom = Fvec2*ones(1,4);

fmat2 = fmat2./denom;

plot(age, fmat2, '-')

hgtmean2 = a.*sum(fmat2,2);

plot(age, hgtmean2, '-', age, hgtmean, '--')

plot(age(2:n), 10.*diff(a.*sum(fmat2,2)))


