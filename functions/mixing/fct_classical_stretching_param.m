function [classical_advection_time1,classical_advection_time2]...
    = fct_classical_stretching_param(model,v)
% Compute the streching times based vorticity and gradient of velocity
%

% Gradient of velocity
grad_w = gradient_mat_2( permute( v , [1 2 4 3]) ,model.grid.dX);
if isfield(model.grid,'mask_keep_cart')
    % Apply mask
    mask_keep_cart = model.grid.mask_keep_cart;
    grad_w = bsxfun(@times,mask_keep_cart,grad_w);
end
% Vorticity
vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);

% Stretching time estimated from vorticity
classical_advection_time1 = 1/std(vort(:));
% Stretching time estimated from vorticity expresse in days
classical_advection_nb_day_1 = classical_advection_time1/3600/24

% Stretching time estimated from velocity gradient
classical_advection_time2 = 1/std(grad_w(:));
% Stretching time estimated from velocity gradient expresse in days
classical_advection_nb_day_2 = classical_advection_time2/3600/24

