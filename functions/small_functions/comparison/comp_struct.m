function b = comp_struct(a,c)

if ~(isstruct(a) && isstruct(c))
    error('They are not structures');
end

fa=fieldnames(a);
la=length(fa);
fc=fieldnames(c);

if la == length(fc);
    b=true;
    for k=1:la
        eval([ 'b = b && comp_field(a.' fa{k} ',c.' fc{k} ');' ]);
    end
else
    b=false;
end

end