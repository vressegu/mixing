%% Test of fct_a

init;

%% Choice of the tracer to be advected
% type_data = 'toy' % Binary tracer filament
 type_data = 'toy2' % smooth tracer filament
% type_data = 'toy3' % smooth tracer blob
% type_data = 'toy4' % smooth tracer filament with varying width
% type_data = 'toy5' % Gaussian random tracer with a spectrum slope -5
% type_data = 'globcurrent' % SST from Globcurrent

% if zoom = true, the spatial domain considered is reduced
% (for Globcurrent data only)
zoom = false;
% if zoom =true, we can better observe specific phenomena
% (e.g. the folding around a vortice)
% but spectra are not accurately computed

%% Choice of the velocity field for the advection

% type_v = 'toy' % rotationally invariant velocity with an angular velocity
%                % which is Gaussian function of the distance at the origin
% type_v = 'toy2' % blob of angular velocity
% type_v = 'toy4' % homogeneous Gaussian tracer with a spectrum slope -6
% type_v = 'toy5' % circular eddy (pure folding)
% type_v = 'toy6' % fat ellispsoidal vortex
% type_v = 'toy7' % Filament of vorticity (pure uniform shear)
 type_v = 'toy8' % ellispsoidal vortex
% type_v = 'toy9' % an ellipsoidal vortex combined with a shear
% type_v = 'globcurrent_v1' % (reinterpolated) geostrophic velocity from
%                           %   v1 SSH data of Globurrent
% type_v = 'globcurrent_v2' % (reinterpolated) geostrophic velocity from
%                           %  v2 SSH data of Globurrent and slighlty
%                           %  filtered

% Note that type_v = toy2,3,4,5,6,7,8 and 9 must be used with one of this
% tracer fields : type_data = toy 2,3,4 or 5
% whereas type_v = globcurrent_v1 and 2 must be used with the tracer
% field type_data = globcurrent

% For the draft paper, we have used the couple
% (type_data = 'toy2' , type_v = 'toy8', 
%       advection_duration_forced.bool = false)
% (type_data = 'globcurrent' , type_v = 'globcurrent_v2' , 
%       advection_duration_forced.bool = false)
% (type_data = 'globcurrent' , type_v = 'globcurrent_v2' , 
%       advection_duration_forced.bool = true)

% Does the velocity field vary in time during the advection?
% (for Globcurrent data only)
v_time_varying = true;
% if v_time_varying = false, only the first time step is used
% if v_time_varying = true, the daily velocity fields are used. Between
% days, the velocity is linearly interpolated in time.

%% Parameters of the advection

% Periodic boundary conditions
periodic_boundary_conditions = false;
% periodic_boundary_conditions is automatically set to true
% for the toy models

% Duration of the advection
advection_duration_forced.bool = true;
%  - If advection_duration_forced.bool = false,
% the duration of the advection is automatically chosen by the algortihm
% in order to make the tracer spectrum slope after the forward advection
% meet a specifc value (see the parameter slope_wanted below).
%  - If advection_duration_forced.bool = true,
% the duration of the advection is specified by the following:
advection_duration_forced.duration =5*24*3600; % 5 days
% Note that if a toy model is used the advection duration is forced

% Value of the wanted tracer spectrum slope after the forward advection
% (if advection_duration_forced.bool = false only)
% slope_wanted = -2.5;
slope_wanted = -3.5;
% tested with satellite images for slope_wanted in the interval [-3.5,2]

% Space step used in the advection
dx =5e3;
% If advection_duration_forced.bool = false,
% the space step dx may be reduced by the algorithm in order to allow for a
% stronger direct cascade of tracer and hence a flatter tracer spectrum
% after the advection

% % dx =5e3;
% % % dx =1e3;
% warning('dx changed');
% dx =1.8e3;

%% Spatial domain and saisonnality
% (if type_data = 'globcurrent' only)


% Suggested study case:

% ACC south of Australia during summer
if zoom % the spatial domain is smaller
    % Boundaries of the spatial domain expressed in longitude/latitude
    lonlat_choose = [ 125 135 -55 -47 ] % ACC south of Australia
%    warning('dx changed');
    %     dx =1e3;
else % the spatial domain is larger
    % Boundaries of the spatial domain expressed in longitude/latitude
    lonlat_choose = [ 105 135 -60 -42 ] % ACC south of Australia
end
% Index of the year in 2011
day =1; % 1st of January 2011 -> Summer in the south hemisphere


% Other tested study cases:

% % ACC south of Australia during winter
% % Boundaries of the spatial domain expressed in longitude/latitude
% lonlat_choose = [ 105 135 -55 -42 ] % ACC south of Australia
% % Index of the year in 2011
% day =201; % 20th of July 2011 -> Winter in the south hemisphere
% % Value of the wanted tracer spectrum slope after the forward advection
% slope_wanted = -3;

% % Gulf stream during winter
% % Boundaries of the spatial domain expressed in longitude/latitude
% lonlat_choose = [ -55 -40 15 39 ]
% % Index of the year in 2011
% day =1; % 1st of January 2011 -> Winter in the north hemisphere
% % Value of the wanted tracer spectrum slope after the forward advection
% slope_wanted = -3;

% % Pacific in the tropics
% % Boundaries of the spatial domain expressed in longitude/latitude
% lonlat_choose = [ -145 -100 -4.5 4.5 ]
% % Index of the year in 2011
% day =1;
% % Value of the wanted tracer spectrum slope after the forward advection
% slope_wanted = -2;


%% Plot paramters

% Regulary plot the spatially filtered tracer during the advection
plot_filtered_field = false;

% If alpha_adv_back=true, the tracer gradient growth-rate is plotted on the
% initial spatial grid x_0 rather than on the final spatial grid x
alpha_adv_back=true;

%% Gather all the parameters in a structure called 'model'
gather_model;

%% Generating the initial tracer field to be advected
[model,sst_ini,sst_masked] = fct_sst_init(model);

% Parameteres of the analytic model for the tracer spectrum tail
model.spectrum_theo0 = fct_param_model_tracer_spectrum_tail(model,sst_ini);

%% Load the the first time step of the velocity field
v = fct_v_init(model);

%% Advection forward
fprintf('Advection Forward \n')

if strcmp(model.type_data(1:min([end 11])),'globcurrent')
    % Sst spectrum before advection
    fct_plot_spectrum_trunc(model,sst_ini,' ini');
end

% Forward advection
[sst_forward,model] = fct_lagrangian_advection4(model, sst_ini,...
    v,sst_masked);

if strcmp(model.type_data(1:min([end 11])),'globcurrent')
    % Sst spectrum after advection
    fct_plot_spectrum_trunc(model,sst_forward,' adv');
    
    % Sst spectrum slope after advection
    slope_sst_forward = fct_estim_spectrum_slope(model,sst_forward)
end

% Save final results
save( [model.folder.folder_simu '/files/Final_all.mat']);

keyboard;