function [alpha2_OW_global,alpha2_OW_local] = ...
    fct_okubo_weiss(model,w,time_adv,day)
% Compute and plot the expected stretching based on the Okubo-Weiss
% assumptions.
% This function can also plot the Okubo-Weiss criterion and
% the ratio r of rotation and strain (when additional_plots = true).
%

% Additional plots ?
additional_plots = false;

if nargin < 4
    day ='0';
end
loc_colorbar = 'southoutside';
colormap_ = 'default';

%% Velocity gradient

% Velocity gradient
grad_w = gradient_mat_2(permute( w,[ 1 2 4 3]),model.grid.dX);
% Mx My 2(grad) 2(w 2D) N_t

% Strain
sigma_s = grad_w(:,:,1,2,:) + grad_w(:,:,2,1,:);
sigma_n = grad_w(:,:,1,1,:) - grad_w(:,:,2,2,:);
vector_sigma_w = sigma_s + 1i * sigma_n;
vector_sigma_w = permute(vector_sigma_w, [ 1 2 3 5 4]); % Mx My 1 N_t
clear sigma_s sigma_n

sigma_w = abs(vector_sigma_w);
sigma_w2 = sigma_w.^2;

% Rotation (vorticity)
vort = grad_w(:,:,1,2,:) - grad_w(:,:,2,1,:);
vort = permute(vort, [ 1 2 3 5 4]);

%% Criterion Q
Q = 1/4 * ( sigma_w2 - vort.^2);

%% Criterion Okubo-Weiss
% ratio of rotation and strain
r_OW = (vort)./sigma_w;

%% Plots

if strcmp(model.type_data,'Gula') && ~model.filtering.smoothing
    caxis([-1 1]*2e-9);
end

% Apply mask
r_OW = model.grid.mask_keep_cart .* r_OW;

% Grid
x = model.grid.x_ref;
y = model.grid.y_ref;

if strcmp(model.type_data(1:min([end 11])),'globcurrent')
    lon = model.grid.lonlat.lon;
    lat = model.grid.lonlat.lat;
    lonlat_ref = model.grid.lonlat.lonlat_ref;
    
    % Inpterpolate to spherical grid
    [LON,LAT]=ndgrid(lon,lat);
    r_OW = ...
        fct_cart2sph(x,y,r_OW,LON,LAT,lonlat_ref);
    sigma_w = ...
        fct_cart2sph(x,y,sigma_w,LON,LAT,lonlat_ref);
    x = lon;y=lat;
end

r2_OW=r_OW.^2;


taille_police = 12;
if strcmp(model.type_data(1:min([end 11])),'globcurrent')
    width=5;
else
    width=4;
end
height=4;
X0=[0 3];

if additional_plots
    figure(48);
    close;
    figure48=figure(48);
    set(figure48,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
    
    imagesc(x,y,(1./r2_OW(:,:,:))');
    colorbar;axis xy;axis equal
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    
    if strcmp(model.type_data(1:min([end 11])),'globcurrent')
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
    else
        ylabel('x',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',taille_police,...
            'FontName','Times')
        xlabel('y)',...
            'interpreter','latex',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'FontName','Times')
    end
    title('OW 1/r2',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    caxis([0 10]);
    %colormap('winter')
    
    colormap(colormap_)
    colorbar('location',loc_colorbar)
    
    % state_OW = 2 : hyperbolic
    % state_OW = 1 : intermediate (stretching linear in time)
    % state_OW = 0 : elliptic
    tol = 1e-1;
    state_OW = 2 * (r2_OW < (1 - tol)) ...
        + 1 * ( ((1-tol) <= r2_OW) & (r2_OW <= (1 + tol)) ) ...
        + 0 * (r2_OW > (1 +tol));
    
    figure(42);
    close;
    figure42=figure(42);
    set(figure42,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
    
    imagesc(x,y,(state_OW(:,:,:))');
    axis xy;axis equal
    caxis([0 2]);
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    if strcmp(model.type_data(1:min([end 11])),'globcurrent')
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
    else
        ylabel('x',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',taille_police,...
            'FontName','Times')
        xlabel('y)',...
            'interpreter','latex',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'FontName','Times')
    end
    title('OW states',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    %colormap('winter')
    
    colormap(colormap_)
    colorbar('location',loc_colorbar)
    
    drawnow
    eval( ['print -depsc ' model.folder.folder_simu '/' ...
        'OW_state.eps']);
end

%% Expected stretching ( = alpha )
% To compute the expected stretching, the compute the expectation of the
% squared tracer gradient norm, assuming that the angle of the intial
% tracer gradient is sampled from random variable uniform in [0, 2*pi]

tol = 1e-1;
stretch = nan(size(r2_OW));

% Elliptic area : |r|>1
iii = (r2_OW > (1 +tol));
r = r_OW(iii);
sqr = sqrt( r.^2 - 1 );
B = sign(r) .* sqr .* sigma_w(iii) * time_adv;
stretch(iii) = cos(B) -1;

% Strain-effective rotation compensated regions : |r|~1
iii = ((1-tol) <= r2_OW) & (r2_OW <= (1 + tol));
r = r_OW(iii);
B = r .* sigma_w(iii) * time_adv;
stretch(iii) = B.^2 * atan(2*pi)/(2*pi) + B * log((2*pi)^2+1)/(2*pi);

% Hyperbolic area : |r|<1
iii = r2_OW < (1 - tol);
r = r_OW(iii);
sqr = sqrt( 1 - r.^2 );
B = sqr .* sigma_w(iii) * time_adv;

% Stretching ( = alpha^2)
stretch(iii) = r .* ( cosh(B)-1 )./(2*pi*sqr) .* ...
    atan( (r-1)./sqr * tanh(2*pi) ) ...
    + 1/(4*pi) * sinh(B) .* log( ( r+cosh(4*pi) ) ./ (r+1) ) ...
    + cosh(B) -1;


figure(3000);
close;
figure3000=figure(3000);
set(figure3000,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');

imagesc(x,y,(stretch)'/time_adv^2);
axis xy;axis equal
% figure;imagesc(x,y,state_OW(:,:,:,day_plot)');axis xy;axis equal
% caxis([0 2]);
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
if strcmp(model.type_data(1:min([end 11])),'globcurrent')
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
else
    ylabel('x',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('y)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
end
title('\hspace{1.5cm}$\alpha_{OW}^2 / t^2$ (Okubo-Weiss)',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')

if strcmp(model.type_data(1:min(end,11)),'globcurrent')
    caxis([0 2e-10]);
end

colormap(colormap_)
colorbar('location',loc_colorbar)

drawnow
eval( ['print -depsc ' model.folder.folder_simu '/' ...
    'OW_stretching' day '.eps']);

alpha2_OW_local = stretch;

% Stretching averaged over space
alpha2_OW_global = mean(stretch(:));

