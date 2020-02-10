function L2 = fct_sq_len_scale_theo(model,spectr_slope,km)
% Compute the square of the length scale
% defined by std(sst)/std(grad(sst))
% assuming that sst as a constant spectrum slope
%

% km = model.grid.km;
% km = 2*pi / min(model.grid.MX.*model.grid.dX);
kinf = pi / max(model.grid.dX);

alpha = - spectr_slope;
kinf = kinf/km;
eta = 1+kinf;

if alpha == 3
    L2 = 1/(2*log(kinf)-3);
elseif alpha == 2
    L2 = 1 / ...
       ( kinf*(eta+1)/eta - 2*log(1+kinf)); 
else
    on_L2 = ...
        2*( 1-eta^(2-alpha) )/(alpha-2) ...
        - kinf*(kinf*(alpha-1)+2)*eta^(1-alpha);
    L2 = (alpha-3) / on_L2;
%     L2 = (alpha-3) / ...
%         ( 2/(alpha-2) ...
%   - (alpha-1) * (kinf)^(-(alpha-3))  );
% %     L2 = (alpha-3)*(alpha-2)/2;
end
L2 = L2 / km^2;