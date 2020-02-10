function f_sph = fct_cart2sph(x,y,f_cart,lon,lat,lonlat_ref,method)
% Reinterpolate data on spherical grid
% The relation between spherical and cartesian coordinates is linearized
%

if nargin < 7
    method = 'cubic';
%     method = 'spline';
end

R = 6367000;
y_sph = R * pi/180*(lat - lonlat_ref(2));
x_sph = R * cos(pi/180*lat) ...
    .* ( pi/180*(lon - lonlat_ref(1)) );
s = size(lon);

f_sph = interp2(x,y,f_cart',x_sph,y_sph,method);
f_sph = reshape(f_sph,s);