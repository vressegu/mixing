function nabla_phi = fct_nabla_phi(model,X)
% This function compute the gradient of the flow
%

X=reshape(X,[model.grid.MX,2]);
% Mx My 2

nabla_phi_x = nan([model.grid.MX,2]);
nabla_phi_y = nan([model.grid.MX,2]);

% Interior
nabla_phi_x(2:end-1,:,:) = 1/(2*model.grid.dX(1)) ...
    *(X(3:end,:,:)-X(1:end-2,:,:));
nabla_phi_y(:,2:end-1,:) = 1/(2*model.grid.dX(2)) ...
    *(X(:,3:end,:)-X(:,1:end-2,:));
% Mx-2 My-2 2 2

% Boundaries
if model.advection.periodic_boundary_conditions
    % Periodic bourndaries
    nabla_phi_x(1,:,:)= 1/(2*model.grid.dX(1)) ...
        *(X(2,:,:)-X(end,:,:));
    nabla_phi_x(end,:,:)= 1/(2*model.grid.dX(1)) ...
        *(X(1,:,:)-X(end-1,:,:));
    nabla_phi_y(:,1,:) = 1/(2*model.grid.dX(2)) ...
        *(X(:,2,:)-X(:,end,:));
    nabla_phi_y(:,end,:) = 1/(2*model.grid.dX(2)) ...
        *(X(:,1,:)-X(:,end-1,:));
else
    % First order approximation
    nabla_phi_x(1,:,:)= 1/model.grid.dX(1) ...
        *(X(2,:,:)-X(1,:,:));
    nabla_phi_x(end,:,:)= 1/model.grid.dX(1) ...
        *(X(end,:,:)-X(end-1,:,:));
    nabla_phi_y(:,1,:) = 1/model.grid.dX(2) ...
        *(X(:,2,:)-X(:,1,:));
    nabla_phi_y(:,end,:) = 1/model.grid.dX(2) ...
        *(X(:,end,:)-X(:,end-1,:));
end

nabla_phi(:,:,:,1) = nabla_phi_x;
nabla_phi(:,:,:,2) = nabla_phi_y;
