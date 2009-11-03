function fdsecobj = section_fd(fdobj, secrng, fdParobj)
basisobj = getbasis(getfd(fdParobj));
nbasis   = getnbasis(basisobj);
fullrng  = getbasisrange(getbasis(fdobj));
if secrng(1) < fullrng(1)-1e-10 || secrng(2) > fullrng(2)+1e-10
    disp([secrng(1)-fullrng(1),secrng(2)-fullrng(2)])
    error('RNG not contained within range of FDOBJ.');
end
n = max(101, 10*nbasis+1);
tvec = linspace(secrng(1),secrng(2),n)';
fmat = eval_fd(tvec, fdobj);
fdsecobj = smooth_basis(tvec, fmat, fdParobj);
fdnames  = getnames(fdobj);
fdsecobj = putnames(fdsecobj, fdnames);