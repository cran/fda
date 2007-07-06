function subbasis = subsref(basisobj, substr)
%  SUBSREF  subscripted reference to a basis object
%    BASIS  ... an object of class 'basis'
%    SUBSTR ... a cell object containing the subscripts

%  last modified 17 December 2007

type = substr.type;

if strcmp(type, '()')
    subs  = substr.subs;
    nsubs = length(subs);
    if nsubs == 1
        dropind = [];
        subs = subs{1};
        for i=1:getnbasis(basisobj);
            if ~any(subs==i)
                dropind = [dropind; i];
            end
        end        
    else
        error('Subscript array is not one-dimensional.');
    end
    subbasis = putdropind(basisobj, dropind);
else
    error('Subscript is not an integer array.');
end

