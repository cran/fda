function rotm = varmx(amat, normalize)
%  VARMX  Does a VARIMAX rotation of a principal components solution
%  Arguments:
%  AMAT      ...  N by K matrix of component loadings
%  NORMALIZE ...  If nonzero, the columns of AMAT are normalized
%                 before computing the rotation matrix.  
%                 The default is 0.
%  Returns:
%  ROTM  ...  Rotation matrix for rotated loadings

%  last modified 27 Oct 2009

if nargin < 2,  normalize = 0;  end

amatd = size(amat);

if length(amatd) ~= 2
    error('AMAT must be two-dimensional')
end

n = amatd(1);
k = amatd(2);
rotm = eye(k);
onek = ones(1,k);
 
if k == 1
    return
end

%  normalize loadings matrix

if normalize
    hvec = sum(amat'.^2)';
    amat = amat./(sqrt(hvec)*onek);
end

eps  = 0.00116;
ccns = 0.7071068;

varold = 0;
varnow = sum(var(amat.^2));

iter = 0;
while abs(varnow - varold) > 1e-7 && iter <= 50
    iter  = iter + 1;
    for j = 1:(k-1)
      for l = (j+1):k
        avecj  = amat(:,j);
        avecl  = amat(:,l);
        uvec   = avecj.^2 - avecl.^2;
        tvec   = 2.*avecj.*avecl;
        aa = sum(uvec);
        bb = sum(tvec);
        cc = sum(uvec.^2 - tvec.^2);
        dd = 2*sum(uvec.*tvec);
        tval = dd - 2.*aa.*bb./n;
        bval = cc - (aa.^2 - bb.^2)./n;

        if tval == bval
          sin4t = ccns;
          cos4t = ccns;
        end

        if tval < bval
          tan4t = abs(tval/bval);
          if tan4t >= eps
            cos4t = 1/sqrt(1+tan4t.^2);
            sin4t = tan4t*cos4t;
          else
            if bval < 0
              sin4t = ccns;
              cos4t = ccns;
            else
              sin4t = 0;
              cos4t = 1;
            end
          end
        end

        if tval > bval
          ctn4t = abs(tval/bval);
          if ctn4t >= eps
            sin4t = 1/sqrt(1+ctn4t.^2);
            cos4t = ctn4t*sin4t;
          else
            sin4t = 1;
            cos4t = 0;
          end
        end

        cos2t = sqrt((1+cos4t)/2);
        sin2t = sin4t/(2*cos2t);
        cost  = sqrt((1+cos2t)/2);
        sint  = sin2t/(2*cost);
        if bval > 0
          cosp = cost;
          sinp = sint;
        else
          cosp = ccns*(cost + sint);
          sinp = ccns*abs(cost - sint);
        end
        if tval <= 0
          sinp = -sinp;
        end

        amat(:,j) =  avecj.*cosp + avecl.*sinp;
        amat(:,l) = -avecj.*sinp + avecl.*cosp;
        rvecj     = rotm(:,j);
        rvecl     = rotm(:,l);
        rotm(:,j) =  rvecj * cosp + rvecl * sinp;
        rotm(:,l) = -rvecj * sinp + rvecl * cosp;

      end
    end

    varold = varnow;
    varnow = sum(var(amat.^2));
end

