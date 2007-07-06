function [fitcell, npar] = bvec2fitcell(bvec, fitcell)
%       transfer parameters from vector BVEC to cell arrays
%                        BWTCELL and AWTCELL

%  Last modified 2 December 2006

fitstruct = fitcell{1};

nvar   = length(fitcell);
Dorder = fitstruct.Dorder;
nbvec  = length(bvec);
npar   = zeros(nvar,Dorder+1);  %  contains numbers of parameters to be estimated.

m2 = 0;
for ivar = 1:nvar
    fitstruct = fitcell{ivar};    
    bwtcell   = fitstruct.bwtcell;   
    awtcell   = fitstruct.awtcell;   
    if isempty(awtcell) || isempty(fitstruct.ufdcell)
        nforce = 0;
    else
        nforce = length(awtcell);
    end
    for jvar=1:nvar
        for ideriv=1:Dorder
            bfdParj = bwtcell{jvar,ideriv};
            if getestimate(bfdParj)
                bfdobjj = getfd(bfdParj);
                bbasisj = getbasis(bfdobjj);
                bnbasis = getnbasis(bbasisj);
                m1 = m2 + 1;
                m2 = m2 + bnbasis;
                if m2 > nbvec
                    error('Number of parameters needed exceeds length of BVEC.');
                end
                npar(ivar,ideriv) = npar(ivar,ideriv) + bnbasis;
                bfdobjj = putcoef(bfdobjj, bvec(m1:m2));
                bwtcell{jvar,ideriv} = putfd(bfdParj, bfdobjj);
            end
        end
    end
    fitstruct.bwtcell = bwtcell;
    
    %  loop through forcing functions to transfer coefficients from
    %  BVEC to the appropriate forcing function weight functions
    
    for k=1:nforce
        afdPark = awtcell{k};
        afdobjk = getfd(afdPark);
        if getestimate(afdPark)
            abasisk = getbasis(afdobjk);
            anbasis = getnbasis(abasisk);
            m1 = m2 + 1;
            m2 = m2 + anbasis;
            if m2 > nbvec
                error('Number of parameters needed exceeds length of BVEC.');
            end
            npar(ivar,Dorder+1) = npar(ivar,Dorder+1) + anbasis;
            afdobjk = putcoef(afdobjk, bvec(m1:m2));
            awtcell{k} = putfd(afdPark, afdobjk);
        end
    end
    fitstruct.awtcell = awtcell;
    fitcell{ivar} = fitstruct;
    
end

if m2 < nbvec
    error('Number of parameters supplied exceeds number required.');
end

