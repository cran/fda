%  test eval_FEM_fd for computing derivatives

p = [0, 0; 1, 0; 0, 1];

t = [1,2,3];

e = [];

order = 2;

basisobj = create_FEM_basis(p, e, t, order);

nfine = 21;
X = linspace(0,1,nfine)';
Y = linspace(0,1,nfine)';

for inode = 1:3*order
    coef = zeros(3*order,1);
    coef(inode) = 1;
    fdobj = fd(coef, basisobj);
    subplot(1,3,1)
    plot_FEM(fdobj, X, Y)
    axis('square')
    title(['Node ',num2str(inode)])
    subplot(1,3,2)
    plot_FEM(fdobj, X, Y, [1,0])
    axis('square')
    title('X-derivative')
    subplot(1,3,3)
    plot_FEM(fdobj, X, Y, [0,1])
    axis('square')
    title('Y-derivative')
    pause
end

if order ==2
for inode = 1:3*order
    coef = zeros(3*order,1);
    coef(inode) = 1;
    fdobj = fd(coef, basisobj);
    subplot(2,2,1)
    plot_FEM(fdobj, X, Y)
    axis('square')
    title(['Node ',num2str(inode)])
    subplot(2,2,2)
    plot_FEM(fdobj, X, Y, [2,0])
    axis('square')
    title('X-second derivative')
    subplot(2,2,3)
    plot_FEM(fdobj, X, Y, [0,2])
    axis('square')
    title('Y-second derivative')
    subplot(2,2,4)
    plot_FEM(fdobj, X, Y, [1,1])
    axis('square')
    title('X-Y-derivative')
    pause
end
end
