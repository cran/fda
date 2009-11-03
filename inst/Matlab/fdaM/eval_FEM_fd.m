function evalmat = eval_FEM_fd(Xvec, Yvec, fdobj, nderivs)
% EVAL_FEM_FD evaluates the FEM fd object at points (Xvec,Yvec)
%
% Arguments:
% Xvec    ... an vector of x-coordinates of length N.  This may also be
%             a matrix with two columns, in which case the first column
%             is assigned to Xvec and the second to Yvec.
% Yvec    ... an vector of y-coordinates of length N.
% FDOBJ   ... A functional data object whose basis is of the FEM type.
% NDERIVS ... A vector of length 2 containing the orders of derivatives
%             for X and Y, respectively.
% 
%        output:
% EVALMAT   an array of the same size as Xvec and Yvec containing the value  
%           of FELSPLOBJ at (Xvec,Yvec).
%
%  Last modified on 27 September 2011 by Jim ramsay.

%  Set up the arguments if the first argument is a matrix with two
%  columns

%  set default values

if nargin < 4,  nderivs = zeros(1,2);  end

if size(Xvec,2) == 2
    if nargin < 2
        error(['First argument is a coordinate matrix and ', ...
               'the second argument is not supplied.']);
    end
    fdobj = Yvec;
    Yvec  = Xvec(:,2);
    Xvec  = Xvec(:,1);
else
    if nargin < 3
        error(['First and second arguments are coordinate vectors ', ...
               'and the third argument is not supplied.']);
    end
end

%  check the type of FDOBJ

if     isa_fd(fdobj)
    if ~strcmp(getbasistype(getbasis(fdobj)), 'FEM')
        error('The basis object for FDOBJ is not of type FEM.');
    end
elseif isa_basis(fdobj)
    if ~strcmp(getbasistype(fdobj), 'FEM')
        error('The basis object for FDOBJ is not of type FEM.');
    end
else
    error('FDOBJ is neither of FD or BASIS class.');
end

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

%  check derivatives

if length(nderivs) ~= 2
    error('NDERIVS not of length 2.');
end
if sum(nderivs) > 2
    error('Maximum derivative order is greater than two.');
end

N = length(Xvec);

%  Augment Xvec and Yvec by one's for computing barycentric coordinates 

Pgpts = [ones(N,1) Xvec Yvec];

%  get basis

if isa_basis(fdobj)
    basisobj = fdobj;
else
    basisobj = getbasis(fdobj);
end

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

if isa_basis(fdobj)
    nbasis  = getnbasis(basisobj);
    coefmat = eye(nbasis);
else
    coefmat = getcoef(fdobj);
end
nsurf   = size(coefmat,2);

% 1st, 2nd, and 3rd vertices of triangles

if     order == 2
    v1 = nodes(nodeindex(:,1),:);
    v2 = nodes(nodeindex(:,3),:);
    v3 = nodes(nodeindex(:,5),:);
elseif order == 1
    v1 = nodes(nodeindex(:,1),:);
    v2 = nodes(nodeindex(:,2),:);
    v3 = nodes(nodeindex(:,3),:);
else
    error('ORDER is neither 1 nor 2.');
end

% denominator of change-of-coordinates change matrix

modJac    = Jvec;
ones3     = ones(3,1);
modJacMat = modJac*ones3';

% 1st, 2nd, and 3rd columns of transformations to barycentric coordinates,
% with a row for each vertex

M1 = [ v2(:,1).*v3(:,2) - v3(:,1).*v2(:,2) ...
       v2(:,2)-v3(:,2)                     ...
       v3(:,1)-v2(:,1)]./modJacMat/2;
M2 = [ v3(:,1).*v1(:,2) - v1(:,1).*v3(:,2) ...
       v3(:,2)-v1(:,2)                     ...
       v1(:,1)-v3(:,1)]./modJacMat/2;
M3 = [ v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2) ...
       v1(:,2)-v2(:,2)                     ...
       v2(:,1)-v1(:,1)]./modJacMat/2;

% Identify triangles containing points in vector (Xvec(i),Yvec(i))
% if no triangle contains a point, ind(i) is NaN

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
                if    sum(nderivs) == 0
                    fval = c1.*(2*baryc1.^2 - baryc1) + ...
                           c2.*(2*baryc2.^2 - baryc2) + ...
                           c3.*(2*baryc3.^2 - baryc3) + ...
                           c4.*(4*baryc1.* baryc2)    + ...
                           c5.*(4*baryc2.* baryc3)    + ...
                           c6.*(4*baryc3.* baryc1);
                elseif nderivs(1) == 1 && nderivs(2) == 0
                    fval = c1.*(4*baryc1 - 1).*M1(indi,2) + ...
                           c2.*(4*baryc2 - 1).*M2(indi,2) + ...
                           c3.*(4*baryc3 - 1).*M3(indi,2) + ...
                           c4.*(4*baryc2.*M1(indi,2)      + ...
                                4*baryc1.*M2(indi,2))     + ...
                           c5.*(4*baryc3.*M2(indi,2)      + ...
                                4*baryc2.*M3(indi,2))     + ...
                           c6.*(4*baryc1.*M3(indi,2)      + ...
                                4*baryc3.*M1(indi,2));
                elseif nderivs(1) == 0 && nderivs(2) == 1
                    fval = c1.*(4*baryc1 - 1).*M1(indi,3) + ...
                           c2.*(4*baryc2 - 1).*M2(indi,3) + ...
                           c3.*(4*baryc3 - 1).*M3(indi,3) + ...
                           c4.*(4*baryc2.*M1(indi,3)      + ...
                                4*baryc1.*M2(indi,3))     + ...
                           c5.*(4*baryc3.*M2(indi,3)      + ...
                                4*baryc2.*M3(indi,3))     + ...
                           c6.*(4*baryc1.*M3(indi,3)      + ...
                                4*baryc3.*M1(indi,3));
                elseif nderivs(1) == 1 && nderivs(2) == 1
                    fval = c1.*(4*M1(indi,2)*M1(indi,3)) + ...
                           c2.*(4*M2(indi,2)*M2(indi,3)) + ...
                           c3.*(4*M3(indi,2)*M3(indi,3)) + ...
                           c4.*(4*M2(indi,2)*M1(indi,3)  + ...
                                4*M2(indi,3)*M1(indi,2)) + ...
                           c5.*(4*M3(indi,2)*M2(indi,3)  + ...
                                4*M3(indi,3)*M2(indi,2)) + ...
                           c6.*(4*M1(indi,2)*M3(indi,3)  + ...
                                4*M1(indi,3)*M3(indi,2));
                elseif nderivs(1) == 2 && nderivs(2) == 0
                    fval = c1.*(4*M1(indi,2)*M1(indi,2)) + ...
                           c2.*(4*M2(indi,2)*M2(indi,2)) + ...
                           c3.*(4*M3(indi,2)*M3(indi,2)) + ...
                           c4.*(8*M2(indi,2)*M1(indi,2)) + ...
                           c5.*(8*M3(indi,2)*M2(indi,2)) + ...
                           c6.*(8*M1(indi,2)*M3(indi,2));
                elseif nderivs(1) == 0 && nderivs(2) == 2
                    fval = c1.*(4*M1(indi,3)*M1(indi,3)) + ...
                           c2.*(4*M2(indi,3)*M2(indi,3)) + ...
                           c3.*(4*M3(indi,3)*M3(indi,3)) + ...
                           c4.*(8*M2(indi,3)*M1(indi,3)) + ...
                           c5.*(8*M3(indi,3)*M2(indi,3)) + ...
                           c6.*(8*M1(indi,3)*M3(indi,3));
                else
                    error('Inadmissible derivative orders.');
                end
                evalmat(i,isurf) = fval;
            else
                c1 = coefmat(nodeindex(indi,1),isurf);
                c2 = coefmat(nodeindex(indi,2),isurf);
                c3 = coefmat(nodeindex(indi,3),isurf);
                if sum(nderivs) == 0
                    fval = (c1.*baryc1 + c2.*baryc2 + c3.*baryc3);
                elseif nderivs(1) == 1 && nderivs(2) == 0
                    fval = c1.*M1(indi,2) + c2.*M2(indi,2) + c3.*M3(indi,2);
                elseif nderivs(1) == 0 && nderivs(2) == 1
                    fval = c1.*M1(indi,3) + c2.*M2(indi,3) + c3.*M3(indi,3);
                else
                    error('Inadmissible derivative orders.');
                end
                evalmat(i,isurf) = fval;
            end
        end
    end
end

