function penaltymat = getbasispenalty (basisobj, Lfdobj)
%GETBASISPENALTY has been replaced by function EVAL_PENALTY.

%  last modified 2 December 2006

warning('Wid:replace', ...
    'GETBASISPENALTY has been replaced by function EVAL_PENALTY.')

%  set up default value for Lfdobj

if nargin < 2 
    Lfdobj = int2Lfd(2);  
end

penaltymat = eval_penalty(basisobj, Lfdobj);