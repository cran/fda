function basismat = getbasismatrix(evalarg, basisobj, nderiv)
%  GETBASISMATRIX   Computes the basis matrix evaluated at arguments in
%    EVALARG associated with basis.fd object BASISOBJ.
%    The returned basis matrix BASISMAT contains the basis
%    derivatives of order NDERIV (0 by default).

%  last modified 11 June 2009 by Jim Ramsay

if nargin < 3,  nderiv = 0;  end

%  check basisobj

if ~isa_basis(basisobj)
    error('Argument BASISOBJ is not a functional basis object');
end

%  search for stored basis matrix

basismat = getbasisvalues(basisobj, evalarg, nderiv);

if ~isempty(basismat)

    %  stored basis matrix found

    return;

else

    %  evaluate basismatrix

    type     = getbasistype(basisobj);
    nbasis   = getnbasis(basisobj);
    params   = getbasispar(basisobj);
    rangeval = getbasisrange(basisobj);
    dropind  = getdropind(basisobj);

    switch type
        case 'bspline'
            rangex   = rangeval;
            if isempty(params)
                breaks = rangex;
            else
                breaks   = [rangex(1), params, rangex(2)];
            end
            norder   = nbasis - length(breaks) + 2;
            basismat = bsplineM(evalarg, breaks, norder, nderiv);
        case 'fourier'
            period   = params(1);
            basismat = fourier(evalarg, nbasis, period, nderiv);
        case 'monom'
            basismat = monomial(evalarg, params, nderiv);
        case 'polyg'
            basismat = polyg(evalarg, params);
        case 'power'
            basismat = powerbasis(evalarg, params, nderiv);
        case 'expon'
            exponents = params;
            basismat = expon(evalarg, exponents, nderiv);
        case 'const'
            basismat = ones(length(evalarg),1);
        case 'QW'
            basismat = QW(evalarg, nderiv);
        case 'QWM'
            basismat = QWM(evalarg, nderiv);
        otherwise
            error('Basis type not recognizable')
    end

    if ~isempty(dropind)
        index = 1:nbasis;
        for i=1:length(dropind)
            index = index(index ~= dropind(i));
        end
        basismat = basismat(:,index);
    end

end