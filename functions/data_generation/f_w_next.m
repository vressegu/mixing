function w_next = f_w_next(model,day_previous,nbp)
% Compute the next veloxity field for the Lagrangian advection method
%

switch model.folder.type_w
    case 'toy'
        w_next = fct_w_toy(model);
    case 'toy2'
        w_next = fct_w_toy2(model); % Pure folding
    case 'toy4'
        w_next = fct_w_toy4(model); % random field
    case 'toy5'
        w_next = fct_w_toy5(model); % Slightly squared vortex
    case 'toy6'
        w_next = fct_w_toy6(model); % Ellipsoid
    case 'toy7'
        w_next = fct_w_toy7(model); % Pure stretching
    case 'toy8'
        w_next = fct_w_toy8(model); % better ellispsoid
    case 'toy9'
        w_next = fct_w_toy9(model); % stretching + ellispsoid
    case 'globcurrent_v1'
        w_next = read_glob_current_w_v1(day_previous ,...
            model.folder.lonlat_choose,model.grid.dX(1));
    case 'globcurrent_v2'
        w_next = read_glob_current_w_v2(day_previous ,...
            model.folder.lonlat_choose,model.grid.dX(1));
    otherwise
        error('unknown type of data');
end
w_next = fct_mirror_on_border(w_next,nbp);
end