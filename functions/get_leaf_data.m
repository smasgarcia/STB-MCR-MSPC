function [data_concat, leaves_used] = get_leaf_data(img_type, date, pot, metadata, Z1_cell, Z2_cell, Z3_cell)
%GET_LEAF_DATA  Pixel spectra of the leaves of one image
%
%   [data_concat, leaves_used] = get_leaf_data(img_type, date, pot, ...
%                                  metadata, Z1_cell, Z2_cell, Z3_cell)
%
%   img_type         : side of the pot, 'G' (left) or 'D' (right)
%   date        : date 'yyyymmdd'
%   pot       : pot number
%   metadata     : struct array from the database (fields fecha = date, tipo = side,
%                  planta = pot, fila/columna = row/column)
%                  locating each image in the cell arrays
%   Z1-3_cell    : cell arrays with the pixel spectra (pixels x bands) of
%                  leaves 1, 2 and 3 of each image
%
%   data_concat : spectra of leaves 1, 2, 3 stacked in this order
%   leaves_used : leaf numbers present (non-empty), in the same order

    idx = find(strcmp({metadata.fecha}, date) & ...
               strcmp({metadata.tipo}, img_type) & ...
               [metadata.planta] == pot);

    data_concat = [];
    leaves_used = [];

    for i = 1:length(idx)
        f = metadata(idx(i)).fila;
        c = metadata(idx(i)).columna;

        Z1 = safe_getfield(Z1_cell, f, c);
        Z2 = safe_getfield(Z2_cell, f, c);
        Z3 = safe_getfield(Z3_cell, f, c);

        if ~isempty(Z1)
            data_concat = [data_concat; Z1];   %#ok<AGROW>
            leaves_used = [leaves_used; 1];    %#ok<AGROW>
        end
        if ~isempty(Z2)
            data_concat = [data_concat; Z2];   %#ok<AGROW>
            leaves_used = [leaves_used; 2];    %#ok<AGROW>
        end
        if ~isempty(Z3)
            data_concat = [data_concat; Z3];   %#ok<AGROW>
            leaves_used = [leaves_used; 3];    %#ok<AGROW>
        end
    end
end
