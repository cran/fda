function penaltymat = getbasispenalty (basisobj, Lfdobj)
%GETBASISPENALTY has been replaced by function EVAL_PENALTY.

%  last modified 14 January 2003

warning('GETBASISPENALTY has been replaced by function EVAL_PENALTY.')

%  set up default value for Lfdobj

if nargin < 2 
    Lfdobj = int2Lfd(2);  
end

penaltymat = eval_penalty(basisobj, Lfdobj);