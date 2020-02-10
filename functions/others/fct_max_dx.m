function  max_dx = fct_max_dx(slope,km)
% Compute the minimal dx to have a slope 'slope'
% and spatial frequency 'km'
%

epsi = 1e-6;

max_dx = pi/km / ( epsi^(1/slope) - 1);
