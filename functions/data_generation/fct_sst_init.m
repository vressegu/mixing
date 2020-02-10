function [model,sst_ini,sst_masked,sst_lonlat] ...
    = fct_sst_init(model)
% Generating the initial tracer field to be advected
%

%% Generating the initial tracer field to be advected
if strcmp(model.type_data(1:min([end 3])),'toy')
    % Spatial grid
    model.grid.MX= 2^8*[1 1];
    model.grid.dX=2e3*[1 1];
    x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
    y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
    model.grid.x=x;
    model.grid.y=y;
    model.grid.x_ref=x;
    model.grid.y_ref=y;
    model.grid.mask_keep_cart = true(model.grid.MX);
    
    % Coriolis frequency
    % (to set up the order of magnitude of the vorticity)
    model.coriolis.f_model='f_plane';
    model.coriolis.earth_radius = 6.4e6;
    angle_grid = 45/360*2*pi; % rad
    OMEGA = 2*pi/(24*60*60); % rad.s^-1
    model.coriolis.f0 = 2* OMEGA * sin( angle_grid );
    
    % For a toy model the advection duration is forced
    advection_duration_forced.bool = true;
end
switch model.type_data
    case 'toy'
        sst_ini = fct_sst_toy(model);
    case 'toy2'
        sst_ini = fct_sst_toy2(model);
    case 'toy3'
        sst_ini = fct_sst_toy3(model);
    case 'toy4'
        sst_ini = fct_sst_toy4(model);
    case 'toy5'
        sst_ini = fct_sst_toy5(model);
    case 'globcurrent'
        % Load SST and interpolate it on a cartesian grid with the
        % prescribed space step dx
        [sst_ini,sst_lonlat,model] = ...
            fct_sst_init_globcurrent(model);
        % This function also generate the spatial grid and 
        % fits parameters for the sst spectrum
        
        if (~model.advection.advection_duration_forced.bool)
            % Test if space step is too large
            % (the space step dx needs to be small enough in order to allow
            %  for a tracer direct cascade stromg enough and hence for
            %  a tracer spectrum flat enough after the advection)
            max_dx = fct_max_dx(...
                model.advection.slope_wanted,model.advection.sst_ini.km);
            if model.grid.dX(1) > max_dx
                % The space step is too large. The sst has to be extracted
                % again and the grid has to be redefined
                
                % Erase the previosu spatial grid
                model=rmfield(model,'grid');
                
                % New space step
                warning(['dx is too large. It is reduced to ' num2str(max_dx)]);
                model.grid.dX = max_dx * ones(1,2);
                
                % Load SST and interpolate it on a cartesian grid with the
                % prescribed space step dx
                [sst_ini,sst_lonlat,model] = ...
                    fct_sst_init_globcurrent(model);
                % This function also generates the spatial grid and
                % fits parameters for the sst spectrum
            end
        end
        clear slope_sst_ini km_sst_ini
        
        % Sst with mask
        sst_masked = sst_ini;
        sst_masked(~model.grid.mask_keep_cart(:)) = -inf;
        sst_masked = reshape(sst_masked,model.grid.MX);
        
        % Make the spherical grid thinner
        model.grid.lonlat.lonref = model.grid.lonlat.lon;
        model.grid.lonlat.latref = model.grid.lonlat.lat;
        lon = model.grid.lonlat.lon;
        lat = model.grid.lonlat.lat;
        lon = linspace(lon(1),lon(end),model.grid.MX(1));
        lat = linspace(lat(1),lat(end),model.grid.MX(2));
        model.grid.lonlat.lon = lon;
        model.grid.lonlat.lat = lat;
        
    otherwise
        error('unkown model.type_data')
end
% Reshape the sst field
sst_ini = reshape( sst_ini ,model.grid.MX);

if strcmp(model.type_data(1:min([end 3])),'toy')
    sst_masked = sst_ini;
end

%%
    function [sst_ini_,sst_lonlat_,model_] = ...
            fct_sst_init_globcurrent(model_)
        % Load SST and interpolate it on a cartesian grid with the
        % prescribed space step model.grid.dX . This function also 
        % generates the spatial grid and fit parameters for the sst spectrum
        %
        
        % Load SST and interpolate it on a cartesian grid with the
        % prescribed space step dx
        [sst_ini_,dx_,lonlat_ref_,mask_keep_cart_, ...
            x_,y_,sst_lonlat_,lon_,lat_] = ...
            read_glob_current_sst(...
            model_.folder.day,model_.folder.lonlat_choose,model_.grid.dX(1));
        
        % Spatial grid
        model_.grid.x_ref = x_;
        model_.grid.y_ref = y_;
        model_.grid.dX = [x_(2)-x_(1) y_(2)-y_(1)];
        model_.grid.x = x_;
        model_.grid.y = y_;
        model_.grid.MX = size(sst_ini_);
        model_.grid.lonlat.lonlat_ref=lonlat_ref_;
        model_.grid.lonlat.lon=lon_;
        model_.grid.lonlat.lat=lat_;
        
        % Mask which delimits the spatial doamin of interest
        model_.grid.mask_keep_cart = mask_keep_cart_;
        model_.grid.mask_ref_ini = mask_keep_cart_;
        
        % Compute the sst spectrum slope and its cutoff frequency
        % (fitting of a paramtetric model for the sst spectrum)
        [slope_sst_ini_,km_sst_ini_] = ...
            fct_estim_spectrum_slope(model_,sst_ini_);
        model_.advection.sst_ini.slope = slope_sst_ini_;
        model_.advection.sst_ini.km = km_sst_ini_;
        
    end

end
