function b = comp_cell(a,c)

if ~(iscell(a) && iscell(c))
    error('They are not cells');
end

la=length(a);

if la == length(c);
    b=true;
    for k=1:la
        b = b && comp_field(a{k},c{k});
    end
else
    b=false;
end


end