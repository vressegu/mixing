function filter_exp = design_filter_5(model)
% Design a Gaussian filter for smoothing between forward and backward
% advection
%

%% Get parameters
sigma_filter=model.advection.sigma_filter;

threshold = 4*sigma_filter;

%% Grid
x= model.grid.x;
x = x -mean(x);
nx = sum(abs(x)<threshold);
y= model.grid.y;
y = y - mean(y);
ny = sum(abs(y)<threshold);
[x,y]=ndgrid(x,y);

iii = (abs(x)<threshold) & (abs(y)<threshold);

%% Filter
filter_exp = exp( - 1/(2*sigma_filter^2) * (x.^2 + y.^2) );

filter_exp = filter_exp(iii);
filter_exp = reshape(filter_exp,[nx ny]);

filter_exp = 1/sum(filter_exp(:)) * filter_exp ;

% figure;imagesc(x(:,1),y(1,:),filter_exp');axis xy;axis equal
