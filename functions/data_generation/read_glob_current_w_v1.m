function [w_cart,dx,lonlat_ref,x,y,w_lonlat,lon,lat] = ...
    read_glob_current_w_v1(day,lonlat_choose,dx)
% Load v1 geostrophic velocity at a specified day of the year 2011
% (default date = 1st of January)
% The SST is reinterpolated on a cartesian grid with a space step dx
% (default dx = 1e4)
%

% Plotting KE and mask in spherical and cartesian coordinates
bool_plot = false;

%% Default parameters
if nargin < 1
    day = 1;
end
if nargin < 2
    lonlat_choose = [ 120 140 -60 -40 ];
end
if nargin <3
    dx = 1e4;
end

%% Date
date_str = fct_data_2011(day);
date_str = [ date_str '000000'];
day_str=sprintf(['%0' num2str(3) 'd'], day);

%% Read velocity globally over the globe for this date

% Folder where data are saved
folder_data = [ pwd '/data/globcurrent/' ...
    'v1_0_global_010_deg_geostrophic_2011/' day_str '/' ...
    date_str '-GLOBCURRENT-L4-CURgeo_0m-ALT_OI-v01.0-fv01.0.nc'];

% Load longitude/latitude grid
lon_global = ncread(folder_data,'lon');
lat_global = ncread(folder_data,'lat');
% Slight modification of the longitude/latitude grid
warning('lon and lat are slightly modified');
dtheta = 0.1;
lon_global = (lon_global(1):dtheta:(lon_global(end)+dtheta))';
lat_global = (lat_global(1):dtheta:(lat_global(end)+dtheta))';


% Load velocity
u = ncread(folder_data,'eastward_geostrophic_current_velocity');
v = ncread(folder_data,'northward_geostrophic_current_velocity');
w_lonlat = u;
w_lonlat(:,:,2) = v;
clear u v

%% Keep only a spatial window
% Center of a sptial window
lonlat_ref = mean(lonlat_choose([[1 2]' [3 4]']),1);
% Earth radius
R = 6367000;
% Mask on the spherical grid for the domain considered
mask_keep_lon = (lonlat_choose(1) <= lon_global) ...
    & (lon_global < lonlat_choose(2));
mask_keep_lat = (lonlat_choose(3) <= lat_global) ...
    & (lat_global < lonlat_choose(4));
mask_keep = bsxfun( @and, mask_keep_lon , mask_keep_lat');

% Spherical coordinates on the domain considered
lon = lon_global(mask_keep_lon);
lat = lat_global(mask_keep_lat);

% Plots
if bool_plot
    MX_keep = [ sum(mask_keep_lon) sum(mask_keep_lat) ];
    for k = 1:2
        w_temp = w_lonlat(:,:,k);
        w_keep(:,:,k) = reshape( w_temp(mask_keep), MX_keep);
    end
    figure;imagesc(lon_global,lat_global,sqrt(sum(w_lonlat.^2,3))');
    axis equal; axis xy;colorbar;
    title('KE the first day')
    figure;imagesc(lon,lat,sqrt(sum(w_keep.^2,3))');
    colorbar;axis equal; axis xy;
    title('masked KE the first day')
end

%% From spherical to cartesian
% Slightly increase the domain for the interpolation on the cartesian grid
lonlat_choose_larger = lonlat_choose + [-4 4 -4 4];

% Compute the domain boundaries in cartesian coordinates
if prod(lonlat_choose_larger(3:4))>=0
    coef = max(cos(pi/180*lonlat_choose_larger(3:4)));
else
    coef = 1;
end
xy_choose = R * [ coef* ...
    pi/180*(lonlat_choose_larger(1:2)-lonlat_ref(1)) ...
    pi/180*(lonlat_choose_larger(3:4)-lonlat_ref(2)) ];

% Cartesian grid
x = xy_choose(1):dx:xy_choose(2);
y = xy_choose(3):dx:xy_choose(4);

% Force the number of points to be even
if mod(length(x),2)== 1
    x(end)=[];
end
if mod(length(y),2)== 1
    y(end)=[];
end
% Update the domain boundaries
xy_choose = [ x(1) x(end) y(1) y(end)];

% Reinterpolate sst on a cartesian grid using a linearized relation
% between spherical and cartesian coordinates
[X,Y]=ndgrid(x,y);
for k = 1:2
    w_cart(:,:,k) = ...
        fct_sph2cart(lon_global,lat_global,w_lonlat(:,:,k),X,Y,lonlat_ref);
end

% Smoothly impose the velocity to zero outside the boundaries
MX_cart = size(X);
mask_smooth = bsxfun( @times, ...
    fct_unity_approx8(MX_cart(1))' , fct_unity_approx8(MX_cart(2)));
w_cart = bsxfun( @times, ...
    mask_smooth , w_cart );

% Check for possible NaNs in the cartesian sst field
if any(isnan(w_cart(:)))
    % Plots cartesian
    figure;imagesc(x,y,any(isnan(w_cart),3)');axis equal; axis xy;
    title('NaNs in the cartesian velocity field')
    error('there is a nan in the cartesian velocity field');
end

% Cartesian plots
if bool_plot
    % Plots cartesian
    figure;imagesc(x,y,sqrt(sum(w_cart.^2,3))');
    colorbar;axis equal; axis xy;
    title('Cartesian KE field')
    figure;imagesc(x,y,mask_smooth');
    colorbar;axis equal; axis xy;
    title('Cartesian mask')
end

%% From cartesian to spherical
if bool_plot
    % Reinterpolate the sst on a spherical grid using a linearized relation
    % between spherical and cartesian coordinates
    [LON,LAT]=ndgrid(lon,lat);
    for k = 1:2
        w_sph(:,:,k) = fct_cart2sph(x,y,w_cart(:,:,k),LON,LAT,lonlat_ref);
    end
    % Plot
    figure;imagesc(lon,lat,(sum(w_sph.^2,3))');axis equal; axis xy;
    colorbar
    title('Spherical KE field')
end
