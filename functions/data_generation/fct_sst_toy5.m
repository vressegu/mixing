function sst=fct_sst_toy5(model)
% Create a Gaussian tracer field with a spectrum slope -5
%

model.odg_b=1;
model.k_inf_on_k1 = 4;
model.slope_b_ini = -5;
sst = init_Spectrum(model);

%% To plot the tracer field
x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
imagesc(x,y,sst');axis xy;
keyboard;