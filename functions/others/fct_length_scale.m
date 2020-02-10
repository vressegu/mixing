function length_scale2 = fct_length_scale(model,f)
% Compute the square of the lenght scale of a scalar field
%

grad_f = gradient_mat_2( f,model.grid.dX);
iii = ( any(isnan(grad_f),3) | any(isinf(abs(grad_f)),3) ...
    | isinf(abs(f)) | ~ model.grid.mask_keep_cart);
iiiv(:,:,1)=iii;
iiiv(:,:,2)=iii;
f_clean = f;
f_clean = f_clean - mean(f_clean(~iii(:)));
f_clean(iii)=0;
f_clean = reshape(f_clean,model.grid.MX);
grad_f_clean = grad_f;
grad_f_clean(iiiv)=0;
grad_f_clean = reshape(grad_f_clean,[model.grid.MX 2]);

norm_T = sum(f_clean(:).^2) .* prod(model.grid.dX);
norm_grad_T = sum(grad_f_clean(:).^2) .* prod(model.grid.dX);

length_scale2 =  norm_T / norm_grad_T ;
