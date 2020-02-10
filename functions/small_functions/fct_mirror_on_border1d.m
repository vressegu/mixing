function f = fct_mirror_on_border1d(f,nbp)
% Recplicate value ouside the domain for nbp pixels
% assuming double periodicity
%

f=[f(end-nbp+1:end,:);f;f(1:nbp,:)];
f=[f(:,end-nbp+1:end) f f(:,1:nbp)];

