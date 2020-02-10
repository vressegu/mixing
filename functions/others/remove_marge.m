function f =remove_marge(f,model)

MX=model.grid.MX;
sf=size(f,2);
% margex=ceil(MX(1)/5);
% margey=ceil(MX(1)/5);
margex=2;
margey=2;
f=reshape(f,[MX sf]);
f=f(1+margex:end-margex,1+margey:end-margey,:);
f=reshape(f,[prod(MX-2*[margex margey]) sf]);