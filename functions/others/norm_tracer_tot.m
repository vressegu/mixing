function [n_grad_T, n_T] = norm_tracer_tot(model, fft_T)
% L^ norm of the tracer and its gradient
%

%% Grid
M=prod(model.grid.MX);
if any( mod(model.grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=model.grid.MX/2;
% if strcmp(meth_anti_alias,'deriv_LowPass')
%     kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ...
%         .* fct_unity_approx5(model.grid.MX(1));
%     ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1]...
%         .* fct_unity_approx5(model.grid.MX(2));
% else
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
% end
[kx,ky]=ndgrid(kx,ky);
k2=(kx/model.grid.dX(1)).^2+(ky/model.grid.dX(2)).^2;

%% Normes

n_grad_T = (2*pi)^2 * k2(:)' * abs(fft_T(:)).^2;
n_T = sum(abs(fft_T(:)).^2);

