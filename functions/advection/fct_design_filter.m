function filter = fct_design_filter(model)
% Design a Gaussian filter which enable to recover the initial smooth field
%

% Get parameters
MX=model.grid.MX;
PX=MX/2;
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end

%% Wave vector
kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
kx=2*pi/model.grid.dX(1)*kx;
ky=2*pi/model.grid.dX(2)*ky;
[kx,ky]=ndgrid(kx,ky);
k=sqrt(kx.^2+ky.^2);
k(PX(1)+1,:)=inf;
k(:,PX(2)+1)=inf;


L2 = model.spectrum_theo0.coef2 - model.spectrum_theo.coef2;
CST = model.spectrum_theo0.coef1 / model.spectrum_theo.coef1;

filter = CST * exp( -1/2 * L2 * k .^2 );
filter = sqrt(filter);
filter(PX(1)+1,:)=0;
filter(:,PX(2)+1)=0;