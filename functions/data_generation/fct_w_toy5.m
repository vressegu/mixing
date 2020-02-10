function w = fct_w_toy5(model)
% Create a circular eddy
%

center1=[0 0];
R1=1;
vort = fct_vortex(model,center1,R1);
w = vort2velocity(model.grid, vort);


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
title('$Vorticity$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'interpreter','latex',...
    'FontName','Times')
axis([x(1) x(end) y(1) y(end)])

% Create the folders
fct_create_folder_plots(model)

eval( ['print -depsc ' model.folder.folder_simu '/vort_w.eps']);


end

function phi = fct_vortex(model_,center,Radius)

x = model_.grid.dX(1)*(0:model_.grid.MX(1)-1);
x=x-mean(x);
x=x-center(1);
y = model_.grid.dX(2)*(0:model_.grid.MX(2)-1);
y=y-mean(y);
y=y-center(2);
r = sqrt(bsxfun(@plus,(x').^2,(y).^2));

sigma = 0.1*sqrt(prod(model_.grid.dX.*model_.grid.MX)); % ~ Rossby radius
radius=0.03*sqrt(prod(model_.grid.dX.*model_.grid.MX));
radius= Radius * radius;

ampli = model_.coriolis.f0/2;
ampli = ampli/8;
theta_dot = ampli * exp(-1/2*(r-radius).^2/sigma^2);

siz = size(theta_dot);
theta_dot = theta_dot(:);
theta_dot(r<radius) = ampli ;
phi = reshape(theta_dot,siz);
end

