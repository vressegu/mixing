function f_cart = fct_sph2cart(lon,lat,f_sph,X,Y,lonlat_ref,method)
% Reinterpolate data on a cartesian grid
% The relation between spherical and cartesian coordinates is linearized
%

if nargin < 7
    method = 'cubic';
%     method = 'spline';
end

R = 6367000;
lat_cart = 1/R * 180/pi*Y + lonlat_ref(2);
lon_cart = 1/R * 180/pi*X ./ cos(pi/180*lat_cart) + lonlat_ref(1);
s = size(lon_cart);

f_cart = interp2(lon,lat,f_sph',lon_cart(:),lat_cart(:),method);
f_cart = reshape(f_cart,s);