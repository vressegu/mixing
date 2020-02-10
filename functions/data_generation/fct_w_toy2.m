function w = fct_w_toy2(model)
% Create an incompressible velocity with a blob of angular velocity 
%

x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
x=x-mean(x);
y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
y=y-mean(y);
r = sqrt(bsxfun(@plus,(x').^2,(y).^2));

sigma = 0.1*sqrt(prod(model.grid.dX.*model.grid.MX)); % ~ Rossby radius
radius=0.03*sqrt(prod(model.grid.dX.*model.grid.MX));
ampli = model.coriolis.f0/2;
ampli = ampli/8;
theta_dot = ampli * exp(-1/2*(r-radius).^2/sigma^2);

siz = size(theta_dot);
r=r(:);
theta_dot = theta_dot(:);
theta_dot(r<radius) = ampli ;
theta_dot = reshape(theta_dot,siz);

[x,y]=ndgrid(x,y);
w(:,:,1) = -y ; w(:,:,2) = x ;
w = bsxfun(@times,w,theta_dot);

%% Plots

vort = vorticity_perso(model.grid, w);

width=3.5;
height=3;
X0 = [0 0];
figure('Name','Angular velocity','NumberTitle','off','Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');
x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
imagesc(x,y,vort');
axis xy; axis equal;
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',11,...
    'FontName','Times')
xlabel('x(m)',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',11,...
    'FontName','Times')
colorbar
title('Vorticity',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'interpreter','latex',...
    'FontName','Times')
axis([x(1) x(end) y(1) y(end)])

% Create the folders
fct_create_folder_plots(model)

eval( ['print -depsc ' model.folder.folder_simu '/vort_w.eps']);

