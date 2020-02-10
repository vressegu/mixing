function fct_plot_vort(model,w)
% Plot the vorticity
%

grad_w = gradient_mat_2( permute( w , [1 2 4 3]) ,model.grid.dX);

if isfield(model.grid,'mask_keep_cart')
    mask_keep_cart = model.grid.mask_keep_cart;
    grad_w = bsxfun(@times,mask_keep_cart,grad_w);
%     w_mask = bsxfun(@times,mask_keep_cart,w);
    w = bsxfun(@times,mask_keep_cart,w);
end

vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);

%% Plot


% Grid
x = model.grid.x_ref;
y = model.grid.y_ref;

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
folder_simu = model.folder.folder_simu;
folder_simu = [ folder_simu '/vort'];
map = 'default';
loc_colorbar = 'southoutside';
% map = model.folder.colormap;

width=3;
% width=9;
% width=12;
height=4;

% width = 13;
% % width = 3.3;
% %     height = 3.2;
ax = [x(end)-x(1) y(end)-y(1)] ;
aspect_ratio = ax(2)/ax(1);
% height = 0.265*aspect_ratio * width;
% % height = 0.26*aspect_ratio * width;
X0=[0 10];

figure(24);
close;
figure24=figure(24);
set(figure24,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');

imagesc(x,y,vort');
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
title('Vorticity ($s^{-1}$)',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
axis([x(1) x(end) y(1) y(end)])
colormap(map)
%colorbar
colorbar('location',loc_colorbar)

hold on;
[X,Y]=ndgrid(x,y);
% [X,Y]=ndgrid(model.grid.x_ref,model.grid.y_ref);
[verts, averts] = streamslice(X',Y',w(:,:,1)',w(:,:,2)','noarrows');
streamline([verts averts]);
% streamslice(X',Y',w(:,:,1)',w(:,:,2)','noarrows');
hold off

drawnow
eval( ['print -depsc ' folder_simu '.eps']);
% keyboard;

%% Plots longitude - latitude


if isfield(model.grid,'lonlat')
    
    
    % Grid
    x = model.grid.x_ref;
    y = model.grid.y_ref;
    lon = model.grid.lonlat.lon;
    lat = model.grid.lonlat.lat;
    lonlat_ref = model.grid.lonlat.lonlat_ref;
    
    [LON,LAT]=ndgrid(lon,lat);
    vort = fct_cart2sph(x,y,vort,LON,LAT,lonlat_ref);
%     w = fct_w_cart2sph(x,y,w,LON,LAT,lonlat_ref);
    
%     lon = lon - mean(lon);
%     lat = lat - mean(lat);
    
%     w(:,:,1) = fct_cart2sph(x,y,w(:,:,1),LON,LAT,lonlat_ref);
%     w(:,:,2) = fct_cart2sph(x,y,w(:,:,2),LON,LAT,lonlat_ref);
    x = lon;y=lat;
    
    % Other parameters
    taille_police = 12;
    %     id_part=1;
    %     type_data = model.type_data;
    folder_simu = model.folder.folder_simu;
    folder_simu = [ folder_simu '/vort_lonlat'];
    map = 'default';
    loc_colorbar = 'southoutside';
    % map = model.folder.colormap;
    
    %         width=9;
    %         % width=12;
    %         height=4;
    
    % width = 13;
    % % width = 3.3;
    % %     height = 3.2;
    ax = [x(end)-x(1) y(end)-y(1)] ;
    aspect_ratio = ax(2)/ax(1);
    % height = 0.265*aspect_ratio * width;
    % % height = 0.26*aspect_ratio * width;
    X0=[0 10];
    
    %% Plot
    
    figure(25);
    close;
    figure25=figure(25);
    set(figure25,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
    
    
    imagesc(x,y,vort');
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel('Lat($^{\circ}$)',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('Lon($^{\circ}$)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    title('Vorticity ($s^{-1}$)',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    axis xy; axis equal
    axis([x(1) x(end) y(1) y(end)])
    colormap(map)
    %colorbar
    colorbar('location',loc_colorbar)
    
    verts_sph = cell(size(verts));
    for k=1:length(verts)
        verts_temp = verts{k};
        R = 6367000;
        verts_temp(:,2) = 1/R * 180/pi*verts_temp(:,2) ...
            + lonlat_ref(2);
        verts_temp(:,1) = ...
            1/R * 180/pi* verts_temp(:,1) ./ cos(pi/180*verts_temp(:,2)) ...
            + lonlat_ref(1);
        verts_sph{k} = verts_temp;
    end
    hold on;
    streamline([verts_sph averts]);
    hold off
    
    drawnow
    
    eval( ['print -depsc ' folder_simu '_lonlat.eps']);
    % keyboard;
    
end

end