function val = safe_getfield(data, arg2, arg3)
%SAFE_GETFIELD  Return a struct field or a cell element, or [] if it does not exist
%
%   val = safe_getfield(S, 'field')   field of struct S ([] if absent)
%   val = safe_getfield(C, row, col)  element C{row, col} ([] if out of range)
    val = [];
    if nargin == 2
        % Struct and field name
        if isstruct(data) && ischar(arg2)
            if isfield(data, arg2)
                val = data.(arg2);
            end
        end
    elseif nargin == 3
        % Cell array and indices
        row = arg2;
        col = arg3;
        if iscell(data) && isscalar(row) && isscalar(col) && ...
           isnumeric(row) && isnumeric(col) && ...
           row > 0 && col > 0 && ...
           size(data,1) >= row && size(data,2) >= col
            val = data{row,col};
        end
    else
        error('safe_getfield requires 2 or 3 input arguments');
    end
end
