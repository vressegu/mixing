function t = fct_unity_approx8(N_t)
% XP must be a 2 x n matrix
% the result is a vector of size n
%

marge = 1;
slop_size_ratio=5;
% slop_size_ratio=20*(N_t/1000);
% slop_size_ratio=20*(N_t/1000);
sslop=ceil(N_t/slop_size_ratio);
N_t=N_t-2*marge;
% N_t=N_t-2*sslop-2*marge;
% N_t=N_t-2*sslop;

t=ones(1,N_t);
t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

t=[zeros(1,marge) t zeros(1,marge)];
% t=[zeros(1,sslop) t zeros(1,sslop)];
