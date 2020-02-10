function [sst_cart,dx,lonlat_ref,mask_keep_cart,x,y,sst_lonlat,lon,lat] = ...
    read_glob_current_sst(day,lonlat_choose,dx)
% Load sst at a specified day of the year 2011
% (default date = 1st of January)
% The SST is reinterpolated on a cartesian grid with a space step dx
% (default dx = 1e4)
%

% Plotting SST in spherical and cartesian coordinates
bool_plot = false;

% Colormap
load('BuYlRd.mat');
map = BuYlRd; clear BuYlRd

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
day_str=sprintf(['%0' num2str(3) 'd'], day);

%% Read SST globally over the globe for this date

% Folder where data are saved
folder_data = [ pwd '/data/globcurrent/' ...
    'third-party_satellite_sea_surface_temperature_l4_odyssea_2011/'];
if day <= 200 % The first and second part of the year are in differnt folders
    folder_data = [folder_data day_str '/'];
else
    folder_data = [folder_data 'next_days/' num2str(day) '_'];
end
folder_data = [ folder_data ...
    date_str '-IFR-L4_GHRSST-SSTfnd-ODYSSEA-GLOB_010-v2.0-fv1.0.nc'];

% Load longitude/latitude grid
lon_global = ncread(folder_data,'lon');
lat_global = ncread(folder_data,'lat');

% Load sst
sst_lonlat = ncread(folder_data,'analysed_sst');

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
    figure;imagesc(lon_global,lat_global,sst_lonlat');
    axis equal; axis xy; colorbar;colormap(map);
    title('Global sst')
    figure;imagesc(lon,lat,reshape(sst_lonlat(mask_keep(:)),...
        [length(lon) length(lat)])');
    axis equal; axis xy;colorbar;colormap(map);
    title('sst to be advected')
    figure;imagesc(lon,lat,mask_keep');axis equal; axis xy;
    title('Mask of the domain');
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
sst_cart = fct_sph2cart(lon_global,lat_global,sst_lonlat,X,Y,lonlat_ref);

% Reinterpolate the mask on a cartesian grid using a linearized relation
% between spherical and cartesian coordinates
mask_keep_ = double(mask_keep);
mask_keep_(~mask_keep)=-inf;
mask_keep_=reshape(mask_keep_,size(sst_lonlat));
mask_keep_cart = ( ...
    fct_sph2cart(lon_global,lat_global,mask_keep_,X,Y,lonlat_ref, ...
    'linear') == 1 );
s = size(mask_keep_cart);
mask_keep_cart(mask_keep_cart == -inf) = false;
mask_keep_cart = reshape(mask_keep_cart,s);

% Check for possible NaNs in the cartesian sst field
if any(isnan(sst_cart(:)))
    % Plots cartesian
    figure;imagesc(x,y,isnan(sst_cart)');axis equal; axis xy;
    title('NaNs in the cartesian sst field')
    error('there is a nan in the cartesian sst field');
end

% Cartesian plots
if bool_plot
    figure;imagesc(x,y,sst_cart');axis equal; axis xy;
    colorbar;colormap(map);
    title('Cartesian sst field')
    figure;imagesc(x,y,mask_keep_cart');axis equal; axis xy;
    title('Cartesian mask')
end

%% From cartesian to spherical
if bool_plot
    % Reinterpolate the sst on a spherical grid using a linearized relation
    % between spherical and cartesian coordinates
    [LON,LAT]=ndgrid(lon,lat);
    sst_sph = fct_cart2sph(x,y,sst_cart,LON,LAT,lonlat_ref);
    % Plot
    figure;imagesc(lon,lat,sst_sph');axis equal; axis xy;
    title('Spherical sst field')
end

