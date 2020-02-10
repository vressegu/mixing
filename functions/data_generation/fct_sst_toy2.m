function sst=fct_sst_toy2(model)
% Create a smooth tracer filament
%

sigma = 0.1;

my2=model.grid.MX(2)/2;
y = 1/2:(my2-1/2);
y  = y / (2*my2);
y = [ -y(end:-1:1) y];

Ly = model.grid.dX(2) * model.grid.MX(2);
y = y * Ly;
sigma = sigma *Ly;

sst = exp( - 1/2 * y.^2 / sigma^2 );

sst=repmat(sst,[model.grid.MX(1) 1]);

%% To plot the tracer field
% x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
% y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
% imagesc(x,y,sst');axis xy;
% keyboard;


