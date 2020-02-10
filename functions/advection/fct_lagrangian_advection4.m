function [T_adv,model,X,X0] = fct_lagrangian_advection4(model, T0, v, varargin)
% Lagrangian advection of T
%

% Plotting KE, vorticity and gradient of velocity
plot_velocity = false;

%% Folder to save plots and files
% Create the folders
fct_create_folder_plots(model)

% Colormap
load('BuYlRd.mat');
model.folder.colormap = BuYlRd; clear BuYlRd

% Create the folders for filtered field
model_filtered = model;
model_filtered.folder.folder_simu = [ model.folder.folder_simu ...
    '/filtered'];
% model_filtered.folder.folder_simu = [ 'images/SQG_MU/' model.type_data '/filtered'];
fct_create_folder_filtered_plots(model_filtered)

%% Time of advection

% Stretching time scale
[model,~,on_tau_local_folding,on_tau_local_stretching,...
    rate_switch] = fct_stretching_time(model,v);

% Classical Eulerian criterion for stretching time (for comparison)
fct_classical_stretching_param(model,v);

% Set the time of advection
if model.advection.advection_duration_forced.bool
    % The advection time is presribed by hand
    model.advection.advection_duration = ...
        model.advection.advection_duration_forced.duration;
    fprintf('The advection time is forced\n');
else
    % Setup the advection time in order to meet a specific spectrum slope
    % (or equivalently a specific correlation length) for the tracer
    model = fct_time_adv_3(model,T0);
    fprintf('The advection time is automatically chosen\n');
end

%% Choice of time step (CFL)
dX=permute(model.grid.dX,[1 3 2]);
bound_dt=sum(bsxfun(@times,abs(v),pi*1./dX),3);
bound_dt=max(bound_dt(:));
bound_dt=1/bound_dt/2;
model.advection.dt_adv = bound_dt;
dt = model.advection.dt_adv;
clear dX dX2 bound_dt

%% Time of advection
N_day = ceil(model.advection.advection_duration/3600/24);
N_t = ceil(N_day*3600*24/dt);
fprintf(['Time of advection : ' num2str(N_t*dt/3600/24) ' days \n']);

%% Grid
% Moving grid
MX=model.grid.MX;
x=model.grid.dX(1)* (0:(MX(1)-1)) ;
y=model.grid.dX(2)* (0:(MX(2)-1)) ;
[x,y]=ndgrid(x,y);
X(:,1)=x(:);
X(:,2)=y(:);
X_forward = X;

%% Preprocessing of the Lagrangian particles, the velocity, the tracer
% and the mask

% Intial velocity field
v_ini = v;

% Create a mask to be advected
mask_adv = double(model.grid.mask_keep_cart);
mask_adv(~model.grid.mask_keep_cart)=-inf;
mask_adv = reshape(mask_adv,MX);

% Number of pixels replicated on the border in order to deal with the
% boundaries, especially in the interpolation procedures
nbp=3;

% Recplicate mask value ouside the domain for nbp pixels
mask_adv = fct_mirror_on_border1d(mask_adv,nbp);
mask_adv0 = mask_adv;

% Recplicate tracer value ouside the domain for nbp pixels
T0 = fct_mirror_on_border1d(T0,nbp);

% Recplicate velocity value ouside the domain for nbp pixels
v = fct_mirror_on_border(v,nbp);

% Reference grid
X0{1}=model.grid.dX(1)* ((-nbp):(MX(1)+nbp-1)) ;
X0{2}=model.grid.dX(2)* ((-nbp):(MX(2)+nbp-1)) ;

% Plotting vorticity, KE and gradient of velocity
if plot_velocity
    fct_plot_vort(model,v_ini);
    fct_plot_KE(model,v_ini);
    fct_plot_grad_v(model,v_ini);
end


%% Initialisation of the advection

% Counts for plots
tt_last=-inf;
vect_t_plot = 0;

% Initialize several vectors

% L2 norm of tracer along time
v_norm_T = model.spectrum_theo0.norm_T0; 
% L2 norm of tracer gradient along time
v_norm_grad_T = model.spectrum_theo0.norm_grad_T0; 
% L2 norm of estimated tracer gradient along time
v_norm_grad_T_estim = model.spectrum_theo0.norm_grad_T0; 
% Mean alpha^2
v_alpha2_m = 0; 
% Mean estimated alpha^2 from our model
v_alpha2_m_estim = 0; 
% Mean estimated alpha^2 from Okubo-Weiss assumptions
v_alpha2_m_estim_OW = 0; 

% Velocity for the first time step
if model.advection.w_time_varying
    
    % Backward advection
    nb_day_adv = N_t*dt/3600/24;
    day_previous = floor(nb_day_adv) + model.folder.day;
    v_previous = f_w_next(model, day_previous, nbp);
    v_next = f_w_next(model, ceil(nb_day_adv) + model.folder.day, nbp);
    t_between_days = nb_day_adv - floor(nb_day_adv);
    
    % Forward advection
    day_previous_For = model.folder.day;
    v_previous_For = v;
    v_next_For = f_w_next(model,day_previous_For+1,nbp);
    clear v
    t_between_days_For = 0;
    
    dt_day = dt/(3600*24);
end

%% Time loop
for t=1:N_t
    
    if model.advection.w_time_varying
        % Estimate v as linear interpolation between two days
        v = (1-t_between_days) * v_previous + t_between_days * v_next;
        % Time relative to the next day
        t_between_days = t_between_days - dt_day;
        
        % Idem for forward advection
        v_For = (1-t_between_days_For) * v_previous_For ...
            + t_between_days_For * v_next_For;
        % Time relative to the last day
        t_between_days_For = t_between_days_For + dt_day;
    end
    
    % Backward advection (to get T(x,t) on a regular grid)
    [X,dX,delta_X_per] = RK4_advection_lagrangienne(model, X, - v, X0);
    
    % Forward advection (to get T(x_0,t) on a regular grid)
    X_forward = RK4_advection_lagrangienne(model,X_forward, v_For, X0);
    
    % Check if there are Nan values
    if any(isnan(X(:)))
        error(['There are NaNs ' ...
            '(possibly because some advected points left the domain)']);
    end
    
    %% PLots and save (each new day)
    tt = floor(t *dt/ (3600*24)); % Number of days
    if tt > tt_last % For each new day
        tt_last = tt;
        
        % Index of the day
        day = num2str(floor(t*dt/(3600*24)));
        
        % Save time index
        vect_t_plot = [vect_t_plot t];
        
        %% Tracer
        
        % Tracer advection
        % (forward advection of the tracer
        % by interpolaton using the backward Lagrangian paths)
        T_adv = interp2(X0{1},X0{2},T0',X(:,1),X(:,2));
        T_adv=reshape(T_adv,MX);
        
        % Mask advection
        mask_adv = interp2(X0{1},X0{2},mask_adv0',X(:,1),X(:,2));
        mask_adv=reshape(mask_adv,MX);
        mask_adv_inf = mask_adv;
        mask_adv(isinf(abs(mask_adv)))=false;
        mask_adv=reshape(mask_adv,MX);
        
        % Check if there are Nan values
        if any(isnan(T_adv(:)))
            error('there are NaNs');
        end
        
        % Plots of the masked tracer
        fct_plot_spatial(model,mask_adv_inf.*T_adv,day);
        if strcmp(model.type_data(1:min([end 11])),'globcurrent')
            % Plots of the masked tracer in spherical coordinates
            fct_plot_spatial_lonlat(model,mask_adv_inf.*T_adv,day);
            
            % Spectrum
            fct_plot_spectrum_trunc(model,mask_adv_inf.*T_adv,day);
            color =['r' 'm' 'c' 'g'];
            % Superimpose the spectrum of the inital spectrum
            for k = 1:length(varargin)
                T_cible = varargin{k};
                fct_plot_spectrum_trunc_cible(model,...
                    T_cible,day,color(k));
            end
        end
        
        %% Statistical description of the tracer
        
        % Compute the parameteres of the analytic model for the tracer
        % spectrum tail
        [model.spectrum_theo,norm_T,norm_grad_T] = ...
            fct_param_model_tracer_spectrum_tail(model,T_adv);
        
        % Estimatation of the L2 norm of the tracer gradient based the
        % inital Eulerian velocity
        norm_grad_T_estim = model.spectrum_theo0.norm_grad_T0 ...
            * ( 1+ (t *dt)^2 * model.spectrum_theo0.on_tau2_global );
        
        % Concatenate these norms at each time in a vector
        v_norm_T = [v_norm_T norm_T];
        v_norm_grad_T = [v_norm_grad_T norm_grad_T];
        v_norm_grad_T_estim = [v_norm_grad_T_estim norm_grad_T_estim];
        
        %% Stretching and Mixing dignostic
        
        % Gradient of the inverse (i.e. backward) flow
        nabla_phi = fct_nabla_phi(model,X + delta_X_per);        
        
        % Mezic criterion
        fct_mezic5(model,nabla_phi,t);
        
        % Backward advection of the flow gradient (for the plots)
        if model.plot.alpha_adv_back
            s=size(nabla_phi(:,:,1,1));
            for i=1:2
                for j=1:2
                    nabla_phi_adv_back=nabla_phi(:,:,i,j);
                    
                    % Recplicate value ouside the domain for nbp pixels
                    nabla_phi_adv_back = ...
                        fct_mirror_on_border1d(nabla_phi_adv_back,nbp);
                    
                    % Backward advection of the flow gradient
                    % (by interpolaton using the forward Lagrangian paths)
                    nabla_phi_adv_back = interp2(X0{1},X0{2},...
                        nabla_phi_adv_back',X_forward(:,1),X_forward(:,2));
                    nabla_phi_adv_back=reshape(nabla_phi_adv_back,s);
                    
                    nabla_phi(:,:,i,j)=nabla_phi_adv_back;
                end
            end
        end
        
        % Apply mask
        if isfield(model.grid,'mask_keep_cart')
            mask_keep_cart = model.grid.mask_keep_cart;
            nabla_phi = bsxfun(@times,mask_keep_cart,nabla_phi);
        end
        
        % Plots carcterisation of stretching in space
        % (alpha^2, beta/alpha, mesochronic vorticity, FTLEs)
        [alpha2_m, cax_alpha] = fct_alpha2(model,nabla_phi,t);
        
        % Estimation of the mean alpha^2 from our model
        alpha2_m_estim = (t *dt)^2 * model.spectrum_theo0.on_tau2_global;
        
        % Okubo Weiss:
        % Estimation of the mean slpha^2 from the Okubo-Weiss assumptions
        % + plots in space
        [alpha2_OW_global,alpha2_OW_local] = fct_okubo_weiss(model,v_ini,...
            t *dt,day);
        
        % Concatenate at each time in vectors for plot along time
        v_alpha2_m = [v_alpha2_m alpha2_m];
        v_alpha2_m_estim = [v_alpha2_m_estim alpha2_m_estim];
        v_alpha2_m_estim_OW = [v_alpha2_m_estim_OW alpha2_OW_global];
        
        
        %% Advection of stretching, folding and shearing time
        
        % Replicate value ouside the domain for nbp pixels
        rate_switch_adv =  fct_mirror_on_border1d(rate_switch,nbp);
        
        % Backward advection of the rate of switching for stretching times
        rate_switch_adv = interp2(X0{1},X0{2},rate_switch_adv',...
            X_forward(:,1),X_forward(:,2));
        rate_switch_adv=reshape(rate_switch_adv,MX);
        
        % Replicate value ouside the domain for nbp pixels
        on_tau_local_folding_adv = ...
            fct_mirror_on_border1d(on_tau_local_folding,nbp);
        
        % Backward advection of the folding time field
        on_tau_local_folding_adv = interp2(X0{1},X0{2},...
            on_tau_local_folding_adv',...
            X_forward(:,1),X_forward(:,2));
        on_tau_local_folding_adv=reshape(on_tau_local_folding_adv,MX);
        
        % Plot the stretching times on a differnt grid
        fct_plot_tau(model, on_tau_local_folding_adv, ...
            on_tau_local_stretching,rate_switch_adv,[],day);
        
        %% Plot of time evolution of tracer and alpha
        fct_plot_tracer_alpha_comp2(model,v_alpha2_m,v_alpha2_m_estim,...
            v_alpha2_m_estim_OW, ...
            v_norm_grad_T/model.spectrum_theo0.norm_grad_T0,...
            v_norm_grad_T_estim/model.spectrum_theo0.norm_grad_T0,...
            dt,vect_t_plot);
        
        
        %% Filtered tracer
        if model.plot.filtered_field
            % Filter isotropic
            filter = fct_design_filter(model);
            fft_T_filtered = fft_T_adv .* filter;
            fct_plot(model_filtered,fft_T_filtered,day);
        end
        
        fprintf([ num2str(t*dt/(24*3600)) ' days of advection \n'])
        
        
        
        %% Save files
        save( [model.folder.folder_simu '/files/' day '.mat'], ...
            'model','t','mask_adv_inf','T_adv','T0','v',...
            'X','X_forward','X0',...
            'v_alpha2_m','v_alpha2_m_estim',...
            'v_norm_grad_T',...
            'v_norm_grad_T_estim');
        
    end
    
    %% Update velocity field
    
    % for backward advection
    if (t_between_days < 0) && model.advection.w_time_varying
        if t < N_t
            day_previous = day_previous - 1;
            v_next = v_previous;
            v_previous = f_w_next(model,day_previous,nbp);
            t_between_days = t_between_days +1;
        end
    end
    
    % for forward advection
    if (t_between_days_For > 1) && model.advection.w_time_varying
        day_previous_For = day_previous_For + 1;
        v_previous_For = v_next_For;
        v_next_For = f_w_next(model,day_previous_For+1,nbp);
        t_between_days_For = t_between_days_For -1;
    end
end

% Tracer interpolation for the forward advection
T_adv = interp2(X0{1},X0{2},T0',X(:,1),X(:,2));
T_adv=reshape(T_adv,MX);

% Mask advection
mask_adv = interp2(X0{1},X0{2},mask_adv0',X(:,1),X(:,2));
mask_adv=reshape(mask_adv,MX);
mask_adv(isinf(abs(mask_adv)))=false;
mask_adv=reshape(mask_adv,MX);
model.grid.mask_keep_cart = mask_adv;

