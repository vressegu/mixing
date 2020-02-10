function X = define_grid(n,BOX)

x = linspace(BOX(1,1),BOX(2,1),n);
y = linspace(BOX(1,2),BOX(2,2),n);
[x,y]=ndgrid(x,y);
