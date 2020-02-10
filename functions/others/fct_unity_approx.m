function obs = fct_unity_approx(obs,x,y,BOX)

s=size(obs);
BOX=BOX';
slop_size_ratio=10;

idxx = (BOX(1,1)<=x & x<=BOX(1,2));
obs(~idxx,:)=0;
idxy = (BOX(2,1)<=y & y<=BOX(2,2));
obs(:,~idxy)=0;

t=ones(sum(idxx),1);
sslop=ceil(sum(idxx)/slop_size_ratio);
t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

obs(idxx,:) = bsxfun(@times, obs(idxx,:), t);

t=ones(1,sum(idxy));
sslop=ceil(sum(idxy)/slop_size_ratio);
t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

obs(:,idxy) = bsxfun(@times, obs(:,idxy), t);
