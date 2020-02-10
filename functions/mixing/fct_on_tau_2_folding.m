function [on_tau2_global, on_tau_local,rate_switch] = ...
    fct_on_tau_2_folding(model,w_eul)
% Compute the inverse of the local folding time
%

% Plot several fields
additional_plot = false;

% Norm of the Eulerian velocity
n_v = sqrt(sum(w_eul.^2,3));

% Found where the velocity is zero
iii = (n_v == 0);
iiiv(:,:,1)=iii;
iiiv(:,:,2)=iii;

if additional_plot
    figure;
    imagesc(model.grid.x_ref,model.grid.y_ref,...
        (model.grid.mask_keep_cart.*n_v)');
    axis xy; axis equal;colormap;
    title(' velocity norm');
end

% Normalized Eulerian velocity
w_eul_n = fct_normalize(w_eul);
% Set the vector to zero when the velocity is zero
s = size(w_eul_n);
w_eul_n(iiiv(:))=0;
w_eul_n = reshape(w_eul_n,s);

% Orthogonal of the normalized Eulerian velocity
w_eul_n_ortho = fct_ortho(w_eul_n);

% Gradient of the normalized Eulerian velocity
grad_w_eul_n = permute( w_eul_n , [1 2 4 3]);
grad_w_eul_n = gradient_mat_2( grad_w_eul_n,model.grid.dX);

% Curvature
grad_w_eul_n = bsxfun(@times, w_eul_n, grad_w_eul_n);
grad_w_eul_n = sum(grad_w_eul_n,3);
grad_w_eul_n = permute( grad_w_eul_n , [1 2 4 3]);
curv = bsxfun(@times, w_eul_n_ortho, grad_w_eul_n);
curv = sum(curv,3);

% Gradient of the Eulerian velocity
grad_w = gradient_mat_2( permute( w_eul , [1 2 4 3]) ,model.grid.dX);

% Applying mask
if isfield(model.grid,'mask_keep_cart')
    mask_keep_cart = model.grid.mask_keep_cart;
    grad_w = bsxfun(@times,mask_keep_cart,grad_w);
    w_eul = bsxfun(@times,mask_keep_cart,w_eul);
end

% Length scale of the velocity to define a bound for the curvature for what
% we call "straight streamlines"
R2 = std(w_eul(:))/(std(grad_w(:)));
R2 = sqrt(3)*pi/2* R2;

if additional_plot
    vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);
    figure;
    imagesc(model.grid.x_ref,model.grid.y_ref,...
        (model.grid.mask_keep_cart.*vort)');
    title('vorticity')
    axis xy; axis equal;colorbar;
    figure;
    imagesc(model.grid.x_ref,model.grid.y_ref,...
        (model.grid.mask_keep_cart.*curv)');
    axis xy; axis equal;colorbar;
    title('curvature')
    figure;
    imagesc(model.grid.x_ref,model.grid.y_ref,...
        (model.grid.mask_keep_cart./curv)');
    axis xy; axis equal;colorbar;
    title('1/curvature')
    caxis([-1 1]*1e5);
end

% Set curvature to zero where the velocity is zero
s = size(curv);
curv(iii)=0;
curv = reshape(curv,s);

% Criterion to discriminate straigth streamlines (where there are mainly
% locally uniform velocity shears) and curved streamlines (where there are 
% mainly folding and shear of angular velocity)
rate_switch = curv * R2;

% Local temporal frequency of the cell
freq = 1/(2*pi) * n_v .* curv;

if additional_plot
    figure;
    imagesc(model.grid.x_ref,model.grid.y_ref,...
        (model.grid.mask_keep_cart.*freq)');
    axis xy; axis equal
    title('frequency and streamlines')
    hold on;
    [X,Y]=ndgrid(model.grid.x_ref,model.grid.y_ref);
    streamslice(X',Y',w_eul(:,:,1)',w_eul(:,:,2)','noarrows');
    hold off
end

% Gradient of the frequency
grad_freq = gradient_mat_2(freq,model.grid.dX);

% Norm of the gradient of the frequency
n_grad_freq = sqrt(sum(grad_freq.^2,3));

if additional_plot
    figure;
    imagesc(model.grid.x_ref,model.grid.y_ref,...
        (model.grid.mask_keep_cart.*n_grad_freq)');
    streamslice(X',Y',w_eul(:,:,1)',w_eul(:,:,2)','noarrows');
    axis xy; axis equal
    caxis([0 5e-11]);
    title('grad(frequency)')
end

% Inverse of the folding time : 1/tau_f
on_tau_local = 1/sqrt(2) * n_v .* n_grad_freq ./ freq;

% Set the inverse of the folding time to zero when the velocity is zero
s = size(on_tau_local);
on_tau_local(iii)=0;
on_tau_local = reshape(on_tau_local,s);

% Only the square of the folding time is used, so we compute its absolute
% value
on_tau_local = abs(on_tau_local);

if isfield(model.grid,'mask_keep_cart')
    % Applying mask and compute the space average of 1/tau_f^2
    on_tau2_global = 1/sum(model.grid.mask_keep_cart(:)) * ...
        sum( model.grid.mask_keep_cart(:) .* on_tau_local(:).^2);
else
    % Compute the space average of 1/tau_f^2
    on_tau2_global = 1/prod(model.grid.MX) * sum(on_tau_local(:).^2);
end
% tau_folding_global = 1/sqrt(on_tau2_global) /(3600*24)


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
    ax = [x(end)-x(1) y(end)-y(1)] ;
    aspect_ratio = ax(2)/ax(1);
    height = aspect_ratio * width;
    X0=[0 10];
    
    figure22=figure(22);
    set(figure22,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
    imagesc(x,y,log((on_tau_local.^2)'));
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
    title('\hspace{0.5cm} $1/\tau_f^2$ ',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    axis xy; axis equal
    colormap(map)
    colorbar
    drawnow
    eval( ['print -depsc ' folder_simu '/on_tau_2_folding.eps']);
end


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