function obs = fct_unity_approx4(XP,BOX)
% XP must be a 2 x n matrix
% the result is a vector of size n
%

XP=double(XP);
x=XP(:,1)';
y=XP(:,2)';
clear XP
obs=zeros(size(x));
BOX=BOX';
% slop_size_ratio=10;
slop_size_ratio=3;
N_t=1000;

idx = (BOX(1,1)<=x & x<=BOX(1,2)) | (BOX(2,1)<=y & y<=BOX(2,2));
% obs(~idx)=0;

t=ones(1,N_t);
sslop=ceil(N_t/slop_size_ratio);
t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

t=[0 t 0];
N_t = N_t+2;
% t=[zeros(1,sslop) t zeros(1,sslop)];
% N_t = N_t+2*sslop;

x_t = linspace(BOX(1,1),BOX(1,2),N_t);
tx=interp1(x_t,t,x);
y_t = linspace(BOX(2,1),BOX(2,2),N_t);
ty=interp1(y_t,t,y);
clear t x_t y_t

obs(idx) =tx.*ty;
