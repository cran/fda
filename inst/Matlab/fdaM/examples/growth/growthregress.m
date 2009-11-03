%  ----------------------------------------------------------
%   Regress velocity for second half against growth for first  
%  ----------------------------------------------------------

%  Last modified 14 November 2008

%  This was an experiment that didn't yield much of interest, but
%  is stored here in case it is revisited.

%  Results for girls:  dividing line is 8 years

n1 = length(age(age <= 7));
indexfine1 = find(agefine <= 7);
hgtbasis1 = create_bspline_basis([1,7], n1+2, 4, age(1:n1));
hgtfdPar1 = fdPar(hgtbasis1, 2, 1e-10);

n2 = length(age(age >= 8));
indexfine2 = find(agefine >=8);
hgtbasis2 = create_bspline_basis([8,18], n2+2, 4, age(age >= 8));
hgtfdPar2 = fdPar(hgtbasis2, 2, 1e-10);

D1hgtfmat1 = eval_fd(agefine(indexfine1), hgtffd, 1);
D1hgtffd1  = smooth_basis(agefine(indexfine1), D1hgtfmat1, hgtfdPar1);
D1hgtffd1  = center(D1hgtffd1);

D2hgtfmat1 = eval_fd(agefine(indexfine1), hgtffd, 2);
D2hgtffd1  = smooth_basis(agefine(indexfine1), D2hgtfmat1, hgtfdPar1);
D2hgtffd1  = center(D2hgtffd1);

D0hgtfmat2 = eval_fd(agefine(indexfine2), hgtffd, 1);
D0hgtffd2  = smooth_basis(agefine(indexfine2), D0hgtfmat2, hgtfdPar2);
D0hgtffd2  = center(D0hgtffd2);

D1hgtfmat2 = eval_fd(agefine(indexfine2), hgtffd, 1);
D1hgtffd2  = smooth_basis(agefine(indexfine2), D1hgtfmat2, hgtfdPar2);
D1hgtffd2  = center(D1hgtffd2);

nbasis1    = 5;
betabasis1 = create_bspline_basis([1,7], nbasis1);
betafdPar1 = fdPar(betabasis1);

nbasis2    = 5;
betabasis2 = create_bspline_basis([8,18], nbasis2);
betafd2    = fd(eye(nbasis2), betabasis2);
betafdPar2 = fdPar(betafd2);

conbasis = create_constant_basis([8,18]);
confdPar = fdPar(conbasis);

betacell = cell(1,nbasis1*nbasis2);
m = 0;
for j=1:nbasis1
    for k=1:nbasis2
        m = m + 1;
        betacell{m} = confdPar;
    end
end

Zfd = D1hgtffd1;
Zstar = inprod(Zfd, betabasis1);

xfdcell  = cell(1,nbasis1*nbasis2);
m = 0;
for j=1:nbasis1
    for k=1:nbasis2
        m = m + 1;
        xfdcell{m}  = fd(Zstar(:,j)', betabasis2(k));
    end
end

Yfd = D0hgtffd2;

fRegressCell = fRegress(Yfd, xfdcell, betacell);

bCell = fRegressCell{4};
Bmat = zeros(nbasis1, nbasis2);
m = 0;
for j=1:nbasis1
    for k=1:nbasis2
        m = m + 1;
        Bmat(j,k) = getcoef(getfd(bCell{m}));
    end
end

betabifdobj = bifd(Bmat, betabasis1, betabasis2);
sfine = linspace(1,7,51)';
tfine = linspace(8,18,51);
betabimat = eval_bifd(sfine, tfine, betabifdobj);

figure(1)
colormap(gray)
surf(sfine, tfine, betabimat')
xlabel('\fontsize{13} Age s')
ylabel('\fontsize{13} Age t')
zlabel('\fontsize{13} \beta(s,t)')
axis([1,7,8,18,-15,20])

figure(2)
contour(sfine, tfine, betabimat')
xlabel('\fontsize{13} Age s')
ylabel('\fontsize{13} Age t')

Yhatfd = fRegressCell{5};

figure(3)
for i=1:ncasef
    plot(Yfd(i));
    lhdl=line(Yhatfd(i));
    set(lhdl, 'color', 'g')
    title(['Case ',num2str(i)])
    axis([8,18,-4,4])
    pause
end

Ymat    = eval_fd(tfine, Yfd);
Yhatmat = eval_fd(tfine, Yhatfd);

Rsqrd = zeros(51,1);
for i=1:51
    SSY = sum(Ymat(i,:).^2);
    SSE = sum((Ymat(i,:)-Yhatmat(i,:)).^2);
    Rsqrd(i) = (SSY-SSE)/SSY;
end

figure(3)
phdl = plot(tfine, Rsqrd, 'b-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Age t')
ylabel('\fontsize{13} R^2(t)')
axis([8,18,0,1])

%  Step through analysis increasing lower limit on covariate age

agelow = 1:4;
Rsqrd = zeros(51,length(agelow));

for age0 = agelow
    disp(['Lower limit for age = ',num2str(age0)])
    index1 = find(age >= age0 & age <= 7);
    n1 = length(age(index1));
    indexfine1 = find(agefine >= age0 & agefine <= 7);
    hgtbasis1  = create_bspline_basis([age0,7], n1+2, 4, age(index1));
    hgtfdPar1  = fdPar(hgtbasis1, 2, 1e-10);
    D1hgtfmat1 = eval_fd(agefine(indexfine1), hgtffd, 1);
    D1hgtffd1  = smooth_basis(agefine(indexfine1), D1hgtfmat1, hgtfdPar1);
    D1hgtffd1  = center(D1hgtffd1);
    nbasis1    = 5;
    betabasis1 = create_bspline_basis([age0,7], nbasis1);
    betafdPar1 = fdPar(betabasis1);
    betacell = cell(1,nbasis1*nbasis2);
    m = 0;
    for j=1:nbasis1
        for k=1:nbasis2
            m = m + 1;
            betacell{m} = confdPar;
        end
    end
    Zfd = D1hgtffd1;
    Zstar = inprod(Zfd, betabasis1);
    xfdcell  = cell(1,nbasis1*nbasis2);
    m = 0;
    for j=1:nbasis1
        for k=1:nbasis2
            m = m + 1;
            xfdcell{m}  = fd(Zstar(:,j)', betabasis2(k));
        end
    end
    fRegressCell = fRegress(Yfd, xfdcell, betacell);
    bCell = fRegressCell{4};
    Bmat = zeros(nbasis1, nbasis2);
    m = 0;
    for j=1:nbasis1
        for k=1:nbasis2
            m = m + 1;
            Bmat(j,k) = getcoef(getfd(bCell{m}));
        end
    end
    betabifdobj = bifd(Bmat, betabasis1, betabasis2);
    sfine = linspace(age0,7,51)';
    tfine = linspace(8,18,51);
    betabimat = eval_bifd(sfine, tfine, betabifdobj);
    figure(1)
    colormap(hot)
    surf(sfine, tfine, betabimat')
    xlabel('\fontsize{13} Age s')
    ylabel('\fontsize{13} Age t')
    zlabel('\fontsize{13} \beta(s,t)')
    axis([1,7,8,18,min(min(betabimat)),max(max(betabimat))])
    figure(2)
    contour(sfine, tfine, betabimat')
    xlabel('\fontsize{13} Age s')
    ylabel('\fontsize{13} Age t')
    Yhatfd = fRegressCell{5};
    Ymat    = eval_fd(tfine, Yfd);
    Yhatmat = eval_fd(tfine, Yhatfd);
    for i=1:51
        SSY = sum(Ymat(i,:).^2);
        SSE = sum((Ymat(i,:)-Yhatmat(i,:)).^2);
        Rsqrd(i,age0) = (SSY-SSE)/SSY;
    end
    pause;
end

figure(3)
phdl = plot(tfine, Rsqrd, 'b-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Age t')
ylabel('\fontsize{13} R^2(t)')
axis([8,18,0,1])

