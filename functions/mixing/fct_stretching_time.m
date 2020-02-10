function [model,on_tau2_global,on_tau_local_folding,on_tau_local_stretching,...
    rate_switch,on_tau_local] = fct_stretching_time(model,v)
% Compute shearing, folding and stretching time
%

% Estimate the inverse of the folding time : 1 / tau_f
[on_tau2_global_folding,on_tau_local_folding,rate_switch] =...
    fct_on_tau_2_folding(model,v);

% Estimate the inverse of the time of uniform shearing : 1 / tau_s
[on_tau2_global_shearing, on_tau_local_stretching] = ...
    fct_on_tau_2_stretching(model,v);

% Applying the mask on the folding and shearing time
if isfield(model.grid,'mask_keep_cart')
    mask_keep_cart = model.grid.mask_keep_cart;
    on_tau_local_folding = mask_keep_cart .* on_tau_local_folding;
    on_tau_local_stretching = mask_keep_cart .* on_tau_local_stretching;
    rate_switch = mask_keep_cart .* rate_switch;
end

% Combine folding and shearing time to obtain a stretching time : 1 / tau
[on_tau2_global, on_tau_local] = fct_plot_tau(model,...
    on_tau_local_folding,on_tau_local_stretching,rate_switch);
model.spectrum_theo0.on_tau2_global=on_tau2_global;
tau_global = 1/sqrt(on_tau2_global)