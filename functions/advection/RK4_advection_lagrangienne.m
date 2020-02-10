function [X,dX,delta_X_per] = RK4_advection_lagrangienne(model, X, w, X0)
% Inegration of the Lagrangian path X with a velcoity w
% using 4th order Runge-Kutta temporal scheme.
% The velocity w is defined on the grid X0
%

% Time step
dt=model.advection.dt_adv;

%% 4th order Runge-Kutta
k1 = deriv_advection_lagrangienne(X, w,X0);
k2 = deriv_advection_lagrangienne( X + k1*dt/2, w,X0);
k3 = deriv_advection_lagrangienne( X + k2*dt/2, w,X0);
k4 = deriv_advection_lagrangienne( X + k3*dt, w,X0);

dX = (dt/3)*(k1/2 + k2 + k3 + k4/2);
X = X + dX;

%% Deal with periodic boundaries conditions
if model.advection.periodic_boundary_conditions
    % Domain size
    Lx=model.grid.dX(1)*model.grid.MX(1);
    Ly=model.grid.dX(2)*model.grid.MX(2);
    
    % Save
    X_true = X;
    
    % Keep X in the domain (in order to enable the evaluation of the
    % velocity at this point in the next iteration)
    X(:,1)= mod(X(:,1),Lx);
    X(:,2)= mod(X(:,2),Ly);
    
    % Save the difference
    % (in order to be able to compute the flow gradient)
    delta_X_per = X_true - X;
    clear X_true
else
    delta_X_per = 0;
end

%% Sub function
    function dX_=deriv_advection_lagrangienne(Xt,wt,X0_)
        % Interpolate the velcoity on the Lagrangian particles
        %
        dX_(:,1)=interp2(X0_{1},X0_{2},wt(:,:,1)',Xt(:,1),Xt(:,2));
        dX_(:,2)=interp2(X0_{1},X0_{2},wt(:,:,2)',Xt(:,1),Xt(:,2));
    end
end