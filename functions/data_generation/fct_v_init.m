function v = fct_v_init(model)
% Load the the first time step of the velocity field
%

if strcmp(model.type_data(1:3),'toy')
    if ~ strcmp(model.folder.type_w(1:3),'toy')
        error(['If the tracer field is a toy model,' ...
            'the velocity field must be a toy model too.']);
    end
elseif strcmp(model.type_data(1:11),'globcurrent')
    if ~ strcmp(model.folder.type_w(1:11),'globcurrent')
        error(['If the tracer field is a satellite image,' ...
            'the velocity field must represent satellite data too.']);
    end
end
switch model.folder.type_w
    case 'toy'
        v = fct_w_toy(model); % rotationally invariant velocity with an angular velocity
        % which is Gaussian function of the distance at the origin
    case 'toy2'
        v = fct_w_toy2(model); % blob of angular velocity
    case 'toy4'
        v = fct_w_toy4(model); % homogeneous Gaussian tracer with a spectrum slope -6
    case 'toy5'
        v = fct_w_toy5(model); % circular eddy (pure folding)
    case 'toy6'
        v = fct_w_toy6(model); % fat ellispsoidal vortex
    case 'toy7'
        v = fct_w_toy7(model); % Filament of vorticity (pure uniform shear)
    case 'toy8'
        v = fct_w_toy8(model); % ellispsoidal vortex
    case 'toy9'
        v = fct_w_toy9(model); % ellipsoidal vortex combined with a shear
    case 'globcurrent_v1'
        v = ...
            read_glob_current_w_v1(model.folder.day,...
            model.folder.lonlat_choose,model.grid.dX(1));
    case 'globcurrent_v2'
        v = ...
            read_glob_current_w_v2_slightFiltered(model.folder.day,...
            model.folder.lonlat_choose,model.grid.dX(1),model);
    otherwise
        error('unknown type of data');
end
