function [Xr,Yr,Zr]=rotate(X,Y,Z)
% X,Y, and Z are the output from the Matlab griddata routine.
% The function ROTATE returns a new set Xr, Yr, and Zr with the property that
% the regression plane they define is parallel to the XY plane.


x=reshape(X,10000,1);y=reshape(Y,10000,1);z=reshape(Z,10000,1);
x=x(~isnan(z));
y=y(~isnan(z));
z=z(~isnan(z));
Xmat=[repmat(1,length(x),1) x y];
pts=(Xmat'*Xmat)\Xmat'*z;
b=pts(2);c=pts(3);
u=cos(atan(-b));v=sin(atan(-b));
Yrot=[u 0 -v;0 1 0;v 0 u];
w=cos(atan(([0 1 0]*Yrot*[b; c; -1])/([0 0 1]*Yrot*[b; c; -1])));
x=sin(atan(([0 1 0]*Yrot*[b; c; -1])/([0 0 1]*Yrot*[b; c; -1])));
Xrot=[1 0 0; 0 w -x; 0 x w];
Rot=Xrot*Yrot;

normu=[reshape(X,10000,1) reshape(Y,10000,1) reshape(Z,10000,1)]*Rot';
Zr=reshape(normu(:,3),100,100);
Yr=reshape(normu(:,2),100,100);
Xr=reshape(normu(:,1),100,100);

%dx=Xr(1,2)-Xr(1,1);dy=Yr(2,1)-Yr(1,1);
%Lapu=del2(Zr,dx,dy);
