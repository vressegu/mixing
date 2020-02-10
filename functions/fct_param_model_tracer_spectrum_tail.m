function [spectrum_theo,norm_T,norm_grad_T] = ...
            fct_param_model_tracer_spectrum_tail(model,T)
% Compute the parameteres of the analytic model for the tracer spectrum 
% tail
%

% Gradient of the intial tracer
grad_T = gradient_mat_2(T,model.grid.dX);

% Deal with nan and inifite values
iii = ( any(isnan(grad_T),3) | any(isinf(abs(grad_T)),3) ...
    | isinf(abs(T)) | ~ model.grid.mask_keep_cart );
iiiv(:,:,1)=iii;
iiiv(:,:,2)=iii;

% Center and remove nan and inf from the initial tracer field
T_clean = T;
T_clean = T_clean - mean(T_clean(~iii));
T_clean(iii)=0;
T_clean = reshape(T_clean,model.grid.MX);

% Remove nan and inf from the initial tracer gradient field
grad_T_clean = grad_T;
grad_T_clean(iiiv)=0;
grad_T_clean = reshape(grad_T_clean,[model.grid.MX 2]);

% L2 norm of initial tracer
norm_T = sum(T_clean(:).^2).* prod(model.grid.dX);
spectrum_theo.norm_T0 = norm_T;

% L2 norm of initial tracer gradient
norm_grad_T = ...
    sum(grad_T_clean(:).^2).*prod(model.grid.dX);
spectrum_theo.norm_grad_T0 = norm_grad_T;

% The two parameters of the analytic model which are expected
% to determine the tracer spectrum tail
spectrum_theo.coef1 = 2 * sqrt( (2*pi)^3 * ...
    spectrum_theo.norm_T0^3 / spectrum_theo.norm_grad_T0 );
spectrum_theo.coef2 =  ...
    spectrum_theo.norm_T0 / spectrum_theo.norm_grad_T0 ;
