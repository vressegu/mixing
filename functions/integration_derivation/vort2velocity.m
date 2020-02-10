function v = vort2velocity(grid, vort)
% Compute the velocity from the vorticity in pseudo spectral
%

%% Fourier transform
fft_f = fft2(vort);

%% Grid
PX=grid.MX/2;
if ~isfield(grid,'k')
    if any( mod(grid.MX,2)~=0)
        error('the number of grid points by axis need to be even');
    end
    kx=1/(grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ...
        .* fct_unity_approx5(grid.MX(1));
    ky=1/(grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1]...
        .* fct_unity_approx5(grid.MX(2));
    [kx,ky]=ndgrid(kx,ky);
    kx=2*pi/grid.dX(1)*kx;
    ky=2*pi/grid.dX(2)*ky;
    k2=kx.^2+ky.^2;
    k2(PX(1)+1,:)=0;
    k2(:,PX(2)+1)=0;
else
    kx=grid.k.kx;
    ky=grid.k.ky;
    k2=grid.k.k2;
end

%% Inversion of the Laplacian
fft_f = bsxfun(@times, - 1./ k2 , fft_f );

%% Gradient
fft_grad_f(:,:,1,:) =  bsxfun(@times, - 1i * ky , fft_f );
fft_grad_f(:,:,2,:) =  bsxfun(@times, + 1i * kx , fft_f);

%% Siungularity and alising
fft_grad_f(1,1,:)=0;
fft_grad_f(PX(1)+1,:,:)=0;
fft_grad_f(:,PX(2)+1,:)=0;

%% Inverse Fourier transform
v = real(ifft2(fft_grad_f));

