function [on_tau2_global, on_tau_local] = ...
    fct_on_tau_2_stretching(model,w_eul)
% Compute the inverse of the local shearing time
%

% Plot several fields
additional_plot = false;

% Norm of the Eulerian velocity
n_v = sqrt(sum(w_eul.^2,3));

% Found where the velocity is zero
iii = (n_v == 0);
iiiv(:,:,1)=iii;
iiiv(:,:,2)=iii;

% Orthogonal of the normalized Eulerian velocity
v_tilde_ortho = fct_ortho(fct_normalize(w_eul));

% Set the orthogonal normalized velocity to zero 
% where the velocity is zero
s = size(v_tilde_ortho);
v_tilde_ortho(iiiv(:))=0;
v_tilde_ortho = reshape(v_tilde_ortho,s);

% Gradient of the norm of the Eulerian velocity
grad_n_v = gradient_mat_2( n_v,model.grid.dX);

% Inverse of the shearing time : 1/tau_s
on_tau_local = 1/sqrt(2) * sum( v_tilde_ortho .* grad_n_v , 3) ;

% Only the square of the shering time is used, so we compute its absolute
% value
on_tau_local = abs(on_tau_local);

if isfield(model.grid,'mask_keep_cart')
    % Applying mask and compute the space average of 1/tau_s^2
    on_tau2_global = 1/sum(model.grid.mask_keep_cart(:)) * ...
        sum( model.grid.mask_keep_cart(:) .* on_tau_local(:).^2);
else
    % Compute the space average of 1/tau_s^2
    on_tau2_global = 1/prod(model.grid.MX) * sum(on_tau_local(:).^2);
end
% tau_stretching_global = 1/sqrt(on_tau2_global) /(3600*24)


%% Plot

if additional_plot
    
    % Grid
    x = model.grid.x_ref;
    y = model.grid.y_ref;
    
    % Other parameters
    taille_police = 12;
    id_part=1;
    type_data = model.type_data;
    folder_simu = model.folder.folder_simu;
    map = 'default';
    % map = model.folder.colormap;
    
    width = 3.3;
    %     height = 3.2;
    ax = [x(end)-x(1) y(end)-y(1)] ;
    aspect_ratio = ax(2)/ax(1);
    height = aspect_ratio * width;
    X0=[0 10];
    
    figure23=figure(23);
    set(figure23,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
    % imagesc(x,y,log((on_tau_local.^2)'));
    imagesc(x,y,((on_tau_local.^2)'));
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel('y(m)',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('x(m)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    title('\hspace{0.5cm} $1/\tau_s^2$ ',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    axis xy; axis equal
    colormap(map)
    colorbar
    drawnow
    eval( ['print -depsc ' folder_simu '/on_tau_2_stretching.eps']);
end

% keyboard;


end

function f_ortho = fct_ortho(f)
% Compute the orthogonal vector in each point of the space
%
f_ortho(:,:,1)= - f(:,:,2);
f_ortho(:,:,2)= + f(:,:,1);
end

function g = fct_normalize(g)
% Compute the orthogonal vector in each point of the space
%
ng = sqrt(sum(g.^2,3));
g = bsxfun( @times, 1./ng , g);
end