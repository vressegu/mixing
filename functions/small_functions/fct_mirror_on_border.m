function w = fct_mirror_on_border(w,nbp)
% Recplicate value ouside the domain for nbp pixels
% assuming zonal periodicity and meridionally constant values on the top
% and bottom
%

w=w(:,:,1)+1i*w(:,:,2);
w=[w(end-nbp+1:end,:);w;w(1:nbp,:)];
w=[repmat(w(:,1),[1 nbp]) w repmat(w(:,end),[1 nbp])];
w(:,:,2)=imag(w);
w(:,:,1)=real(w(:,:,1));
