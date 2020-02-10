function w = fct_w_toy9(model)
% Create an ellipsoidal vortex combined with a shear
%

% Vortex
c_g = 0*[1 1];
R=1.5;
excentricity = 1.7;
R1=R;
vort = fct_vortex(model,c_g,R1,excentricity);

% Shear
c_g = 0;
R=3;
coef = 0.1;
center1=c_g*[1 1];
R1=R;
vort_shear = coef*fct_shear(model,center1,R1);

% Total vorticity
w = vort2velocity(model.grid, vort + vort_shear );

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

function vort = fct_vortex(model_,center,Radius,excentricity)

LX=model_.grid.MX .* model_.grid.dX;

vort = fct_vortex_(model_,center,Radius,excentricity);
vort = vort + fct_vortex_(model_,center+2*[1 0],Radius,excentricity);
vort = vort + fct_vortex_(model_,center-2*[1 0],Radius,excentricity);
vort = vort + fct_vortex_(model_,center+2*[1 1],Radius,excentricity);
vort = vort + fct_vortex_(model_,center-2*[1 1],Radius,excentricity);
vort = vort + fct_vortex_(model_,center+2*[1 -1],Radius,excentricity);
vort = vort + fct_vortex_(model_,center+2*[-1 +1],Radius,excentricity);
vort = vort + fct_vortex_(model_,center+2*[0 1],Radius,excentricity);
vort = vort + fct_vortex_(model_,center-2*[0 1],Radius,excentricity);

end

function phi = fct_vortex_(model_,center,Radius,excentricity)

LX = model_.grid.dX .* model_.grid.MX;
x = model_.grid.dX(1)*(0:model_.grid.MX(1)-1);
x=x-mean(x);
x=x-LX(1)/2*center(1);
y = model_.grid.dX(2)*(0:model_.grid.MX(2)-1);
y=y-mean(y);
y=y-LX(2)/2*center(2);
[x,y]=ndgrid(x,y);
r = sqrt((excentricity*(x+y)).^2 + (y-x).^2 );

sigma = 0.1*sqrt(prod(model_.grid.dX.*model_.grid.MX)); % ~ Rossby radius
sigma = Radius * sigma;
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


function vort = fct_shear(model_,center,Radius)

LX=model_.grid.MX .* model_.grid.dX;

vort = fct_shear_(model_,center,Radius);
vort = vort + fct_shear_(model_,center+2*[1 0],Radius);
vort = vort + fct_shear_(model_,center-2*[1 0],Radius);
vort = vort + fct_shear_(model_,center+2*[1 1],Radius);
vort = vort + fct_shear_(model_,center-2*[1 1],Radius);
vort = vort + fct_shear_(model_,center+2*[1 -1],Radius);
vort = vort + fct_shear_(model_,center+2*[-1 +1],Radius);
vort = vort + fct_shear_(model_,center+2*[0 1],Radius);
vort = vort + fct_shear_(model_,center-2*[0 1],Radius);

end

function phi = fct_shear_(model_,center,Radius)

LX = model_.grid.dX .* model_.grid.MX;
x = model_.grid.dX(1)*(0:model_.grid.MX(1)-1);
x=x-mean(x);
x=x-LX(1)/2*center(1);
y = model_.grid.dX(2)*(0:model_.grid.MX(2)-1);
y=y-mean(y);
y=y-LX(2)/2*center(2);
[x,y]=ndgrid(x,y);
r=abs(x);

sigma = 0.1*sqrt(prod(model_.grid.dX.*model_.grid.MX)); % ~ Rossby radius
sigma= Radius * sigma;

ampli = model_.coriolis.f0/2;
ampli = ampli/8;
theta_dot = ampli * exp(-1/2*r.^2/sigma^2);

siz = size(theta_dot);
theta_dot = theta_dot(:);
phi = reshape(theta_dot,siz);
end


