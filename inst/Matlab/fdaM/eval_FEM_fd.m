function evalmat = eval_FEM_fd(Xvec,Yvec,fdobj)
% EVAL_FEM_FD evaluates the FEM fd object at points (Xvec,Yvec)
%
%        arguments:
% FELSPLOBJ a FELspline object
% Xvec         an array of x-coordinates.
% Yvec         an array of y-coordinates.
% SOLUTION  a vector containing a solution to the smoothing problem defined
%           in FELSPLOBJ.  Must be the same length as FELSPLOBJ.SOLUTION.
% 
%        output:
% EVALMAT   an array of the same size as Xvec and Yvec containing the value of 
%           FELSPLOBJ at (Xvec,Yvec).

%  Last modified on 21 August 2010 by Jim ramsay.

%  check Xvec

if ~isa(Xvec,'double')
   error('Xvec is not a numerical array')
else
   Xvec  = Xvec(:);     % treat Xvec as a column vector
end

%  check Yvec

if ~isa(Yvec,'double')
   error('Yvec is not a numerical array')
elseif length(Yvec(:))~=length(Xvec)
   error('Yvec is not the same length as Xvec')
else
   Yvec=Yvec(:);     % treat Yvec as a column vector
end

N = length(Xvec);

%  matrix of local coordinates 

Pgpts = [ones(N,1) Xvec Yvec];

%  get basis

basisobj = getbasis(fdobj);

%  get nodes and index

params    = getbasispar(basisobj);
p         = params.p;
t         = params.t;
t         = t(:,1:3);
order     = params.order;
nodes     = params.nodes;
nodeindex = params.nodeindex;
Jvec      = params.J;

%  get the coefficient matrix

coefmat = getcoef(fdobj);
nsurf   = size(coefmat,2);

% 1st, 2nd, and 3rd vertices of triangles

if     order == 2
    pg1 = nodes(nodeindex(:,1),:);
    pg2 = nodes(nodeindex(:,3),:);
    pg3 = nodes(nodeindex(:,5),:);
elseif order == 1
    pg1 = nodes(nodeindex(:,1),:);
    pg2 = nodes(nodeindex(:,2),:);
    pg3 = nodes(nodeindex(:,3),:);
else
    error('ORDER is neither 1 nor 2.');
end

% denominator of change-of-coordinates change matrix

modJac = Jvec;
ones3 = ones(3,1);
modJacMat = modJac*ones3';

% 1st, 2nd, and 3rd rows of change-of-coordinates matrix

M1 = [(pg2(:,1).*pg3(:,2))-(pg3(:,1).*pg2(:,2)) ...
       pg2(:,2)-pg3(:,2)                        ...
       pg3(:,1)-pg2(:,1)]./modJacMat;
M2 = [(pg3(:,1).*pg1(:,2))-(pg1(:,1).*pg3(:,2)) ...
       pg3(:,2)-pg1(:,2)                        ...
       pg1(:,1)-pg3(:,1)]./modJacMat;
M3 = [(pg1(:,1).*pg2(:,2))-(pg2(:,1).*pg1(:,2)) ...
       pg1(:,2)-pg2(:,2)                        ...
       pg2(:,1)-pg1(:,1)]./modJacMat;

% identify element containing point in vector (Xvec(i),Yvec(i))
% if no element contains a point, ind(i) is NaN

tricoef = tricoefCal(p, t);
ind = zeros(N,1);
for i=1:N
   ind(i) = insideIndex(Xvec(i), Yvec(i), p, t, tricoef);
end

%  interpolate values

evalmat = NaN.*zeros(N, nsurf);

for isurf=1:nsurf
    for i=1:N
        indi = ind(i);
        if ~isnan(indi)
            %  change to barycentric coordinates
            baryc1 = (M1(indi,:).*Pgpts(i,:))*ones3;
            baryc2 = (M2(indi,:).*Pgpts(i,:))*ones3;
            baryc3 = (M3(indi,:).*Pgpts(i,:))*ones3;
            if order == 2
                c1 = coefmat(nodeindex(indi,1),isurf);
                c2 = coefmat(nodeindex(indi,3),isurf);
                c3 = coefmat(nodeindex(indi,5),isurf);
                c4 = coefmat(nodeindex(indi,2),isurf);
                c5 = coefmat(nodeindex(indi,4),isurf);
                c6 = coefmat(nodeindex(indi,6),isurf);
                fval = c1.*(2*baryc1.^2 - baryc1) + ...
                       c2.*(2*baryc2.^2 - baryc2) + ...
                       c3.*(2*baryc3.^2 - baryc3) + ...
                       c4.*(4*baryc1.* baryc2) + ...
                       c5.*(4*baryc2.* baryc3) + ...
                       c6.*(4*baryc3.* baryc1);
                evalmat(i,isurf) = fval;
            else
                c1 = coefmat(nodeindex(indi,1),isurf);
                c2 = coefmat(nodeindex(indi,2),isurf);
                c3 = coefmat(nodeindex(indi,3),isurf);
                fval = c1.*baryc1 + c2.*baryc2 + c3.*baryc3;
                evalmat(i,isurf) = fval;
            end
        end
    end
end

