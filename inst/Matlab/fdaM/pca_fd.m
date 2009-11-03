function pcastr = pca_fd(fdobj, nharm, harmfdPar, centerfns)
%  PCA Functional principal components analysis with regularization
%
%  Arguments:
%  FDOBJ     ... Functional data object (a struct object)
%  NHARM     ... Number of principal components to be kept. Default 2
%  HARMFDPAR ... A functional parameter object specifying the
%                basis, differential operator, and level of smoothing
%                for eigenfunctions.
%  CENTERFNS ... If 1, the mean function is first subtracted from each 
%                function.  1 is the default.
%
%  Returns:
%  A struct object PCASTR with the fields:
%  HARMFD  ... A functional data object for the harmonics or eigenfunctions
%  VALUES  ... The complete set of eigenvalues
%  HARMSCR ... A matrix of scores on the principal components or harmonics
%  VARPROP ... A vector giving the proportion of variance explained
%                 by each eigenfunction
%  FDHATFD ... A functional data object for the approximation to the
%              FDOBJ based on NHARM principal components
%  MEANFD  ... A functional data object giving the mean function
%
%  If NHARM = 0, all fields except MEANFD are empty.

%  Last modified:  24 April 2012 by Jim Ramsay

%  check FDOBJ

if ~isa_fd(fdobj)
    error ('First argument is not a functional data object.');
end

%  get basis information for functional data

fdbasis  = getbasis(fdobj);

%  set up default values

if nargin < 4
    centerfns = 1;   %  subtract mean from data before PCA
end

if nargin < 3
    %  default Lfd object: penalize 2nd deriv., lambda = 0
    Lfdobj    = int2Lfd(2);
    lambda    = 0;
    harmfdPar = fdPar(fdbasis, Lfdobj, lambda);
else
    %  check harmfdPar object
    if ~isa_fdPar(harmfdPar)
        if isa_fd(harmfdPar) || isa_basis(harmfdPar)
            harmfdPar = fdPar(harmfdPar);
        else
            error(['HARMFDPAR is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
        end
    end
end

if nargin < 2
    nharm = 2;  %  default to two harmonics
end

%  compute mean function

meanfd = mean(fdobj);

if nharm == 0
    pcastr.harmfd  = [];
    pcastr.values  = [];
    pcastr.harmscr = [];
    pcastr.varprop = [];
    pcastr.fdhatfd = [];
    pcastr.meanfd  = meanfd;
    return
end

%  --------------------   begin principal components analysis  ------------

if centerfns ~= 0
    % center data
    fdobj = center(fdobj);
end

%  set up HARMBASIS

harmbasis = getbasis(getfd(harmfdPar));
nhbasis = getnbasis(harmbasis);

%  set up LFDOBJ

Lfdobj = getLfd(harmfdPar);
Lfdobj = int2Lfd(Lfdobj);

%  set up LAMBDA

lambda = getlambda(harmfdPar);

%  get coefficient matrix for FDOBJ and its dimensions

coef   = getcoef(fdobj);
coefd  = size(coef);
nbasis = coefd(1);
nrep   = coefd(2);
ndim   = length(coefd);

if nrep < 2
    error('PCA not possible without replications');
end

%  compute CTEMP whose cross product is needed

if ndim == 3
    nvar  = coefd(3);
    ctemp = zeros(nvar*nbasis,nrep);
    for j = 1:nvar
        index = (1:nbasis) + (j-1)*nbasis;
        ctemp(index,:) = coef(:,:,j);
    end
else
    nvar = 1;
    ctemp = coef;
end

%  set up cross product Lmat for harmonic basis, 
%  roughness penalty matrix Rmat, and
%  penalized cross product matrix Lmat.

Lmat = eval_penalty(harmbasis, int2Lfd(0));
if lambda > 0
    Rmat = eval_penalty(harmbasis, Lfdobj);
    Lmat = Lmat + lambda .* Rmat;
end
Lmat = (Lmat + Lmat')/2;

%  Choleski factor Mmat of Lmat = Mmat'*Mmat

Mmat    = chol(Lmat);
Mmatinv = inv(Mmat);

%  coefficient cross product matrix for covariance operator

Wmat = ctemp*ctemp'./nrep;

%  set up matrix for eigenanalysis depending on whether
%  a special basis was supplied for the eigenfunctions or not

Jmat = inprod(harmbasis, fdbasis);
MIJW = Mmatinv'*Jmat;

if nvar == 1
    Cmat = MIJW*Wmat*MIJW';
else
    Cmat = zeros(nvar*nhbasis);
    for i = 1:nvar
        indexi =   (1:nbasis) + (i-1)*nbasis;
        for j = 1:nvar
            indexj = (1:nbasis) + (j-1)*nbasis;
            Cmat(indexi,indexj) = MIJW*Wmat(indexi,indexj)*MIJW';
        end
    end
end

% Eigenanalysis

Cmat = (Cmat + Cmat')./2;
[eigvecs, eigvals] = eig(Cmat);
[eigvals, indsrt ] = sort(diag(eigvals));
eigvecs = eigvecs(:,indsrt);
neig    = nvar*getnbasis(harmbasis);
indx    = neig + 1 - (1:nharm);
eigvals = eigvals(neig + 1 - (1:neig));
eigvecs = eigvecs(:,indx);
sumvecs = sum(eigvecs);
eigvecs(:,sumvecs < 0) = -eigvecs(:,sumvecs < 0);

varprop = eigvals(1:nharm)./sum(eigvals);

%  Set up fdnames for harmfd

harmnames = getnames(fdobj);
%  Name and labels for harmonics
harmlabels = ['I   '; 'II  '; 'III '; 'IV  '; 'V   '; ...
              'VI  '; 'VII '; 'VIII'; 'IX  '; 'X   '];
if nharm <= 10
    harmnames2    = cell(1,2);
    harmnames2{1} = 'Harmonics';
    harmnames2{2} = harmlabels(1:nharm,:);
    harmnames{2}  = harmnames2;
else
    harmnames{2} = 'Harmonics';
end
%  Name and labels for variables
if iscell(harmnames{3})
    harmnames3    = harmnames{3};
    harmnames3{1} = ['Harmonics for ',harmnames3{1}];
    harmnames{3}  = harmnames3;
else
    if ischar(harmnames{3}) && size(harmnames{3},1) == 1
        harmnames{3} = ['Harmonics for ',harmnames{3}];
    else
        harmnames{3} = 'Harmonics';
    end
end

%  set up harmfd

if nvar == 1
    harmcoef = Lmat\eigvecs;
else
    harmcoef = zeros(nbasis,nharm,nvar);
    for j = 1:nvar
        index     = (1:nbasis) + (j-1)*nbasis;
        eigvecsj  = eigvecs(index,:);
        harmcoef(:,:,j) = Lmat\eigvecsj;
    end
end

harmfd = fd(harmcoef, harmbasis, harmnames);
% harmfd = fd(harmcoef, fdbasis);

%  set up harmscr

if nvar == 1
    harmscr = inprod(fdobj, harmfd);
else
    harmscr       = zeros(nrep, nharm, nvar);
    coefarray     = getcoef(fdobj);
    harmcoefarray = getcoef(harmfd);
    for j=1:nvar
        coefj     = squeeze(coefarray(:,:,j));
        harmcoefj = squeeze(harmcoefarray(:,:,j));
        fdobjj    = fd(coefj, fdbasis);
        harmfdj   = fd(harmcoefj, fdbasis);
        harmscr(:,:,j) = inprod(fdobjj,harmfdj);
    end
end

%  set up functional data object for fit to the data

if nvar == 1
    fdhatcoef = harmcoef*harmscr';
else
    fdhatcoef = zeros(nbasis,nrep,nvar);
    for j=1:nvar
        fdhatcoef(:,:,j) = harmcoef(:,:,j)*harmscr(:,:,j)';
    end
end
fdhatfd = fd(fdhatcoef, harmbasis);

%  set up structure object PCASTR

pcastr.harmfd  = harmfd;
pcastr.values  = eigvals;
pcastr.harmscr = harmscr;
pcastr.varprop = varprop;
pcastr.fdhatfd = fdhatfd;
pcastr.meanfd  = meanfd;

