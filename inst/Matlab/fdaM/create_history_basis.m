function bibasisobj = create_history_basis(M, B, tn, DomainType)
%CREATE_HISTORY_BASIS computes indices of the three nodes 
%   defining each triangular element covering either a triangle
%   or parallelogram domain for a set of points (s,t) having the
%   property t - \delta <= s <= t.  If the domain is triangular, that
%   the lower limit is the maximum of 0 and t - \delta.  This domain
%   is used for the historical linear model that approximates the 
%   value y(t) of a functional dependent variable on the basis of
%   \int x(s,t)\beta(s,t) ds where s is restricted in this way.
%   The basis system for approximating bivariate regression function
%   \beta is a set of piecewise linear "hat" or "tent" functions
%   with hexagonal support constructed from these triangular elements.
%   It is assumed that y(t) is observed for times 0 <= t <= TN.

%   The number of basis functions in this system is equal to the 
%   number of vertices defining the triangular elements, called the
%   NODES of the system.  This is the number of rows of the output
%   array ELENODES that provides the indices of the triangles with
%   vertices at each node.  
%
% Inputs :
%
% M  ... number of intervals into which the domain 
%         boundaries are divided.
% TN ... time of the last observation of the Y's.
% B  ... width of the domain of integration in terms
%         of number elements.  Default value is M.
% DomainType ... Type of domain as follows
%      DomainType = 1:  Triangular,    domain of s = [ 0,T]
%  In the following figures, elements are lettered and nodes are numbered.
%  Both elements and nodes are numbered row-wise from left to right
%    and from the bottom to the top.
%
%          4--5--6
%          |b/|d/
%          |/c|/ 
%          2--3
%          |a/
%          |/
%          1
%
%      DomainType = 2:  Parallelogram, domain of s = [-T,T]
%  For M = 2 and B = 2, there are 8 elements and 9 nodes.
%
%          7--8--9
%         /|f/|h/
%        /e|/g|/ 
%       4--5--6
%      /|b/|d/
%     /a|/c|/
%    1--2--3
%
%  The total number of elements and nodes are:
%      DomainType = 1:  Elements ... 2*B(M-1), and nodes ... (B+1)(M-1+B/2).
%      DomainType = 2:  Elements ... 2*M*B,    and nodes ... (M+1)*(B+1).
%
%  However, when B = 0 and lag = 0, the two domain types are equivalent, 
%    the number of elements is M and the number of nodes is M+1. 
%  The default value for DomainType is 1.
%
% Outputs :
% BASISOBJ:  a functional data object with
%    TYPE   ... 'history'
%    RANGE  ... [0,TN]
%    NBASIS ... number of nodes
%    PARAMS ... An array whose first NBASIS*3 elements are the 
%               array eleNodes described below  stored column-wise, 
%               and whose remaining elements (an even number)
%               are, first, the equally spaced arguments s defining
%               the mesh, and second, the same number of 
%               equally spaced arguments t.  The number of these
%               arguments depends on the geometry of the domain
%               and on the number of basis functions in a somewhat
%               complicated way.
%    
% eleNodes:  a matrix with 
%          number rows = total number of elements
%          number cols = 3
% Each row contains the node numbers defining the element, numbered starting
%    with the upper left node and going anti-clockwise.
%
% However, when B = 0, the first column contains the left node, 
%                      the second the right node, and the third column is 0.
%
% For the above DomainType = 2 diagram, eleNodes returns
%   
%  4  1  2
%  4  2  5
%  5  2  3
%  5  3  6
%  7  4  5
%  7  5  8
%  8  5  6
%  8  6  9

%  Last modified:  20 July 2006

% LAG : difference between M and the index of the upper element.
%  Nonzero lags are not allowed in this function.  
%       See functions NodeIndexation and ParalleloGrid for
%       possibly positive lages.

lag = 0;  

%  set some default parameter values

if nargin < 4, DomainType = 1; end

if DomainType ~= 1 && DomainType ~= 2
    error('Argument DomainType has incorrect value');
end

if nargin < 3, lag = 0; end
lag = floor(lag);
if lag < 0 || lag >= M
    error('Lag is not between 0 and M - 1');
end

if nargin < 2, B = M; end

if B < 0 || B > M
    error('Argument B has incorrect value');
end
    
Bp1 = B + 1;

rangeval = [0, tn];

type = 'history';

%  Domain [tn, tn]:

if B == 0
    nbasis = M;
    eleNodes = zeros(nbasis,3);
    eleNodes(:,1) = (1:(M)  )';
    eleNodes(:,2) = (2:(M+1))';
    Si = linspace(0,tn,M+1)';
    Ti = Si;
    params = [reshape(eleNodes,nbasis*3,1); Si; Ti];
    bibasisobj = bibasis(type, rangeval, rangeval, nbasis, params);
else
    
    %  Domain [0, tn]
    
    if DomainType == 1
        nbasis    = M^2 - (M-B)^2;
        eleNodes = zeros(nbasis,3);
        ielem = 0;
        mm3   = 0;
        mm4   = 1;
        for m = 1:B
            %  first element in row is right
            ielem = ielem + 1;
            eleNodes(ielem,:) = [mm4+1,mm3+1,mm4+2];
            for b = 1:(m-1)
                bp1 = b + 1;
                %  left element
                ielem = ielem + 1;
                eleNodes(ielem,:) = [mm4+bp1,mm3+b,  mm3+bp1];
                %  right element
                ielem = ielem + 1;
                eleNodes(ielem,:) = [mm4+bp1,mm3+bp1,mm4+b+2];
            end
            mm3 = mm3 + m;
            mm4 = mm4 + m + 1;
        end
        for m = (B+1):M
            for b = 1:B
                bp1 = b + 1;
                %  left element
                ielem = ielem + 1;
                eleNodes(ielem,:) = [mm4+b,mm3+b,  mm3+bp1];
                %  right element
                ielem = ielem + 1;
                eleNodes(ielem,:) = [mm4+b,mm3+bp1,mm4+bp1];
            end
            mm3 = mm3 + B + 1;
            mm4 = mm4 + B + 1;
        end
        Si = 0;
        Ti = 0;
        for m=1:B
            Si = [Si, linspace(0, m*lambda, m+1)];
            Ti = [Ti, (m*lambda).*ones(1,m+1)];
        end
        for m = (B+1):M
            Si = [Si, linspace((m-B)*lambda, m*lambda, B+1)];
            Ti = [Ti, (m*lambda).*ones(1,B+1)];
        end
        params = [reshape(eleNodes,nbasis*3,1); Si; Ti];
        bibasisobj = bibasis(type, rangeval, rangeval, nbasis, params);
    end
    
    %  Domain [-tn,tn]:
    
    if DomainType == 2
        nbasis = 2*M*B;
        eleNodes = zeros(nbasis,3);
        for m = 1:M
            mm1 = m - 1;
            m2  = 2*mm1*B;
            m3  = mm1*Bp1;
            m4  = m*Bp1;
            for b = 1:B
                bp1 = b + 1;
                %  left element
                b2  = 2*b + m2 - 1;
                eleNodes(b2,:) = [m4+b,m3+b  ,m3+bp1];
                %  right element
                b2  = b2 + 1;
                eleNodes(b2,:) = [m4+b,m3+bp1,m4+bp1];
            end
        end
        s0 = (-lag - B)*lambda;  %  beginning of interval of integration
        is = linspace(s0,-lag*lambda,B+1);  %  abscissa boundaries of elements
        Si = is;
        Ti = zeros(1,B+1);
        for m = 1:M
            Si = [Si,is+m*lambda];
            Ti = [Ti,ones(1,B+1)*m*lambda];
        end
        params = [reshape(eleNodes,nbasis*3,1); Si; Ti];
        bibasisobj = bibasis(type, rangeval, rangeval, nbasis, params);
    end
    
end

