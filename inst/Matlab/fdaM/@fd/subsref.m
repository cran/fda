function subfd = subsref(fd, substr)
%  SUBSREF  subscripted reference to a functional data object
%    FD     ... an object of class 'fd'
%    SUBSTR ... a cell object containing the subscripts

%  last modified 20 July 2006

coef = getcoef(fd);

type = substr.type;

if strcmp(type, '.')
    subfd = eval(['fd',substr.subs{2}]);
elseif strcmp(type, '()')
    sizec = size(coef);
    ndim  = length(sizec);
    subs  = substr.subs;
    nsubs = length(subs);
    switch nsubs
        case 1
            switch ndim
                case 1
                    error('Too many subscripts.');
                case 2
                    newcoef = coef(:,subs{1});
                case 3
                    newcoef = coef(:,subs{1},:);
                otherwise
                    error('Wrong no. dimensions for COEF.');
            end
        case 2
            switch ndim
                case 1
                    error('Too many subscripts.');
                case 2
                    error('Too many subscripts.');
                case 3
                    newcoef = coef(:,subs{1},subs{2});
                otherwise
                    error('Wrong no. dimensions for COEF.');
            end
        case 3
            error('Too many subscripts.');
    end

    subfd.coef     = newcoef;
    subfd.basisobj = getbasis(fd);
    subfd.fdnames  = getnames(fd);

    subfd = class(subfd, 'fd');
else
    error('Cell subscripting is not allowed, refer by field name.');
end

