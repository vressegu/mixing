%% Gather all the parameters in a structure called 'model'
lonlat_choose_str = fct_vect2str(lonlat_choose);
day_str = fct_data_2011(day);
model.folder.day = day;
model.grid.dX=dx;
model.folder.lonlat_choose=lonlat_choose;
model.folder.main_folder_simu = [ 'images/T_' type_data ...
    '_w_' type_v ...
    '_lonlat_' lonlat_choose_str ...
    '_day_' day_str ...
    '_wanted_slope_' num2str(slope_wanted)];
model.plot.alpha_adv_back = alpha_adv_back;
% model.mirror=mirror;clear mirror
model.plot.filtered_field=plot_filtered_field;clear plot_filtered_field
if ~advection_duration_forced.bool
    model.advection.wanted_scale = ['slope' num2str(slope_wanted)];
    model.advection.slope_wanted = slope_wanted;
    % model.advection.wanted_scale = wanted_scale;
end
% model.centered= centered;
model.advection.advection_duration_forced=advection_duration_forced;
model.physical_constant.f0 = 1e-4;
model.folder.folder_simu = [ 'images/Lagrangian_advection/T_' type_data '_v_' type_v ];
% model.advection.sigma_filter=sigma_filter;
model.type_data=type_data;
model.folder.type_w=type_v;
model.advection.w_time_varying = v_time_varying;
if strcmp(model.type_data(1:min([end 3])),'toy')
    periodic_boundary_conditions = true;
end
model.advection.periodic_boundary_conditions = periodic_boundary_conditions;

% Colormap
load('BuYlRd.mat');
model.folder.colormap = BuYlRd; clear BuYlRd