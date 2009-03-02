function basismat = getbasismatrix(evalarg, basisobj, nderiv)
%  GETBASISMATRIX   Computes the basis matrix evaluated at arguments in
%    EVALARG associated with basis.fd object BASISOBJ.
%    The returned basis matrix BASISMAT contains the basis
%    derivatives of order NDERIV (0 by default).

%  last modified 1 November 2007

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
    nbasis   = nbasis + length(dropind);

    switch type
        case 'bspline'
            rangex   = rangeval;
            if isempty(params)
                basismat = monomial(evalarg, 0:nbasis-1, nderiv);
            else
                breaks   = [rangex(1), params, rangex(2)];
                norder   = nbasis - length(breaks) + 2;
                basismat = bsplineM(evalarg, breaks, norder, nderiv);
            end
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
        otherwise
            error('Basis type not recognizable')
    end

    if length(dropind) > 0
        index = 1:nbasis;
        for i=1:length(dropind)
            index = index(index ~= dropind(i));
        end
        basismat = basismat(:,index);
    end

end