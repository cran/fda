function Lfdobj = vec2Lfd(bwtvec, rangeval)
%  VEC2LFD converts a vector of length m to a linear differential operator 
%  object of order m.  The range of the functional data object in any cell 
%  is set to RANGEVAL (defaults to [0 1]) and a constant basis is used.
%
%  In the event that BWTVEC is already a linear differential operator
%  object, it returns the object.  
%
%  Returns a Lfd object defining a constant-coefficient linear 
%  differential operator. 

%  Last modified 15 May, 2005. 

if  isa_Lfd(bwtvec),
    Lfdobj = bwtvec;
else
    if(nargin<2), rangeval = [0 1]; end
    %  check BWTVEC

    if ~isnumeric(bwtvec)
        error('Argument not a vector and not a linear differential operator.');
    end
    m = length(bwtvec);

    %  set up the list object for the homogeneous part

    if m==0
        %  if derivative is zero, BWTLIST is empty
        bwtlist = [];
    else
        conbasis = create_constant_basis(rangeval);
        bwtcell  = cell(m,1);
        for j = 1:m,
            bwtcell{j} = fdPar(fd(bwtvec(j), conbasis));
        end
    end

    Lfdobj = Lfd(m, bwtcell);
end
