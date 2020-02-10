function sst=fct_sst_toy4(model)
% Create a smooth tracer filament with varying width
%

sigma = 0.1;

my2=model.grid.MX(2)/2;
y = 1/2:(my2-1/2);
y  = y / (2*my2);
y = [ -y(end:-1:1) y];

mx2=model.grid.MX(1)/2;
x = 1/2:(mx2-1/2);
x  = x / (2*mx2);
x = [ -x(end:-1:1) x];

[x,y]=ndgrid(x,y);

Ly = model.grid.dX(2) * model.grid.MX(2);
y = y * Ly;
sigma = sigma *Ly;

sigmax = 0.1;
sigma = sigma * (1 -1/3 * exp(-1/2 * (x-0.2).^2/sigmax^2) );

sst = exp( - 1/2 * y.^2 ./ (sigma.^2) );

%% To plot the tracer field
% x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
% y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
% imagesc(x,y,sst');axis xy;
% keyboard;
