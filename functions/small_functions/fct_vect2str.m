function vect_str = fct_vect2str(v)
% Create a string from the vector v, adding '_' between each coefficient of
% the vector
%

vect_str = [];
for k=1:length(v)
    vect_str = [ vect_str '_' num2str(v(k))];
end
vect_str(1)=[];