function fct_plot_spatial(model,T_adv_part,day)
% This function creates some plot online and save it
%

%% Get paramters

% Grid
x = model.grid.x_ref;
y = model.grid.y_ref;
% lon = model.grid.lonlat.lon;
% lat = model.grid.lonlat.lat;
% lonlat_ref = model.grid.lonlat.lonlat_ref;
% 
% [LON,LAT]=ndgrid(lon,lat);
% T_adv_part_lonlat = ...
%     fct_cart2sph(x,y,T_adv_part,LON,LAT,lonlat_ref);

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
folder_simu = model.folder.folder_simu;
plot_moments = false;
map = model.folder.colormap;

%% One particle
X0=[0 0];

if ( (eval(day) == 0) && ...
        strcmp(model.type_data,'Perturbed_vortices') )
    width = 3.2;
    height = 3.2;
    figure1=figure(1);
    set(figure1,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
    contourf(x,y,T_adv_part');
    x= model.grid.dX(1)*(0:model.grid.MX(1)-1);
    y= model.grid.dX(2)*(0:model.grid.MX(2)-1);
    Lx = model.grid.dX(1) * model.grid.MX(1);
    sigma= 2 * Lx/15;
    center1x=x(1/4*model.grid.MX(1)+1);
    center1y=y(1/4*model.grid.MX(2)+1);
    nsig=40;
    dist = 1.5;
    rate = 0.3;
    sigma = sigma/nsig;
    center1x= 2e4;
    center1y=y(1/4*model.grid.MX(2)+1);
    coord1=[center1x center1y];
    size_square = 10e4;
    redline1 = [ [coord1(1)-size_square/2 coord1(2)-size_square/2] ; ...
        [coord1(1)+size_square/2 coord1(2)-size_square/2] ; ...
        [coord1(1)+size_square/2 coord1(2)+size_square/2] ; ...
        [coord1(1)-size_square/2 coord1(2)+size_square/2] ; ...
        [coord1(1)-size_square/2 coord1(2)-size_square/2] ];
    hold on;
    if strcmp(model.type_data,'Perturbed_vortices')
        plot(redline1(:,1),redline1(:,2),'r','LineWidth',3);
    elseif strcmp(model.type_data,'spot6')
        plot(redline1(:,1),redline1(:,2),'b','LineWidth',3);
    else
        error('wrong type of data?');
    end
    center1x= x(1/2*model.grid.MX(1)+1) - 2e4;
    center1y=y(3/4*model.grid.MX(2)+1);
    coord1=[center1x center1y];
    redline1 = [ [coord1(1)-size_square/2 coord1(2)-size_square/2] ; ...
        [coord1(1)+size_square/2 coord1(2)-size_square/2] ; ...
        [coord1(1)+size_square/2 coord1(2)+size_square/2] ; ...
        [coord1(1)-size_square/2 coord1(2)+size_square/2] ; ...
        [coord1(1)-size_square/2 coord1(2)-size_square/2] ];
    plot(redline1(:,1),redline1(:,2),'b','LineWidth',3);
    hold off;
    
else
    width = 3.3;
    
%     height = 3.2;
    
    ax = [x(end)-x(1) y(end)-y(1)] ;
    aspect_ratio = ax(2)/ax(1);
    height = aspect_ratio * width;
    
    figure1=figure(1);
    set(figure1,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
%     imagesc(lon,lat,T_adv_part_lonlat');
    imagesc(x,y,T_adv_part');
    
end

% caxis([-1 1]*1e-3);
% caxis([0 1]);
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
title(['\hspace{0.5cm} $t=' num2str(day) '$ day '],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
colormap(map)
colorbar
drawnow
eval( ['print -depsc ' folder_simu '/one_realization/'...
    num2str(day) '.eps']);       


end

