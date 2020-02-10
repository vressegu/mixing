function [model,time_adv] = fct_time_adv_3(model,sst_ini)
% Estimate the advection time to obtain a specific tracer spectrum slope
% ( model.advection.slope_wanted) after the advection
%

% Theoretical sst squared correlation length after advection
% assuming a specific parametric form for the sst spectrum and
% assuming that the cutoff frequency km_sst_ini (planetary scale)
% remains constant during the advection
L2_wanted = fct_sq_len_scale_theo(model,...
    model.advection.slope_wanted ,model.advection.sst_ini.km);

% Required (theoretical) sst correlation length after advection
length_scale_wanted = sqrt(L2_wanted);

% square of the initial sst correlation lenght
L2_sst_init = fct_length_scale(model,sst_ini);

% Width of the Gaussian filter that could be used to go
% from the required sst correlation lenght after advection
% to the initial sst correlation lenght
sigma_filter = sqrt( L2_sst_init - length_scale_wanted^2);
model.advection.sigma_filter = sigma_filter;
% This quantity is used to specify the time of advection

% Estimate the needed (time-dependent) growth rate alpha^2
alpha2_wanted = model.spectrum_theo0.coef2 / (sigma_filter^2) - 1 ;
alpha2_wanted = max([0 alpha2_wanted]);

% Estimate the advection time
time_adv2 = 1/ (alpha2_wanted * model.spectrum_theo0.on_tau2_global );
time_adv = sqrt(time_adv2);
model.advection.advection_duration = time_adv;

