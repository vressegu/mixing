function sst=fct_sst_toy3(model)
% Create a smooth tracer blob
%

odg_b = 1e-3;

mx2=model.grid.MX(1)/2;
my2=model.grid.MX(2)/2;
sigma = 0.3;
y = 1/2:(my2-1/2);
y  = y /my2;
y = [ -y(end:-1:1) y];
x = 1/2:(mx2-1/2);
x  = x /mx2;
x = [ -x(end:-1:1) x];
[x,y]=ndgrid(x,y);
sst = exp( - 1/2 * (x.^2 + y.^2) / sigma^2 );

sst = odg_b * (2*sst-1);

%% To plot the tracer field

% x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
% y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
% imagesc(x,y,sst');axis xy;
% hold on;
% keyboard;