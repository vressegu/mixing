function res = norm_table(X,p)
if nargin == 1
    p=2;
end
d= length(size(X));
X=abs(X).^p;
for k=1:d
    X = sum(X);
end
res=X^(1/p);
end