
function b =comp_field(d,e)
if isstruct(d)
    if isstruct(e)
        b=comp_struct(d,e);
    else
        b=false;
    end
elseif islogical(d)
    if islogical(e)
        d=d(:);
        e=e(:);
        if length(d)==length(e)
            b=all(eq(d,e));
        else
            b=false;
        end
    else
        b=false;
    end
elseif isnumeric(d)
    if isnumeric(e)
        d=d(:);
        e=e(:);
        if length(d)==length(e)
            b=all(eq(d,e));
        else
            b=false;
        end
    else
        b=false;
    end
elseif isa(d,'stk_dataframe')
    if isa(e,'stk_dataframe')
        b=isequal(d,e);
%         b=all(b(:));
%         d=d(:);
%         e=e(:);
%         if length(d)==length(e)
%             b=all(isequal(d,e));
%         else
%             b=false;
%         end
    else
        b=false;
    end
elseif ischar(d)
    if ischar(e)
        b=strcmp(d,e);
    else
        b=false;
    end
elseif iscell(d)
    if iscell(e)
        d=d(:);
        e=e(:);
        b=comp_cell(d,e);
    else
        b=false;
    end
else
    error('not coded yet for this type of data');
end

end