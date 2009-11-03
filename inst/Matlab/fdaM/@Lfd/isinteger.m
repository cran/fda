function intwrd = isinteger(Lfdobj)
% ISINTEGER returns 1 of WFD and AFD are both zero functions.

%  Last modified 3 January 2008

%  check WFDCELL for emptyness or all zero

wfdcell = getwfdcell(Lfdobj);
wintwrd = 1;
if ~isempty(wfdcell)
    nderiv = Lfdobj.nderiv;
    for j=1:nderiv
        wfdParj = wfdcell{j};
        wfdj    = getfd(wfdParj);
        if any(getcoef(wfdj) ~= 0.0)
            wintwrd = 0;
        end
    end
end

intwrd = wintwrd;