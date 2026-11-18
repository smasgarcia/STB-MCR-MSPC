function img = rebuild_map_image(pot, date, side, ...
    metadataG, metadataD, ...
    G_pixin1, G_pixin2, G_pixin3, G_pixout, ...
    D_pixin1, D_pixin2, D_pixin3, D_pixout, ...
    info_table, T2_limit)
%REBUILD_MAP_IMAGE  Pixel map of one image: green = below the limit, red = above
%
%   img = rebuild_map_image(pot, date, side, metadataG, metadataD, ...
%         G_pixin1, G_pixin2, G_pixin3, G_pixout, ...
%         D_pixin1, D_pixin2, D_pixin3, D_pixout, info_table, T2_limit)
%
%   pot     : pot number
%   date      : date 'yyyymmdd'
%   side       : 'G' (left) or 'D' (right) side of the pot
%   metadataG/D: struct arrays from the database (fields fecha = date, planta = pot,
%                fila/columna = row/column) that
%                locate each image in the pixel-index cell arrays
%   *_pixin1-3 : cell arrays with the linear indices of the pixels of leaves
%                1-3 in the resized image (548 x 313 pixels)
%   info_table : table with columns pot, date, leaf and T2_all (per-pixel
%                T2 of each leaf, in the same order as the pixel indices)
%   T2_limit   : limit (1 when T2_all is already divided by the limit, as in
%                the table produced by s05)
%
%   Returns an RGB image (548 x 313 x 3): leaf pixels with T2 <= limit in
%   green and pixels with T2 > limit in red; background in black.

    % size of the resized images (rows x columns); the pixel indices
    % (pixin) are linear indices in an image of this size
    nRows = 313;
    nCols = 548;

    if strcmpi(side, 'G')
        meta = metadataG;
        pixin = {G_pixin1, G_pixin2, G_pixin3};
    else
        meta = metadataD;
        pixin = {D_pixin1, D_pixin2, D_pixin3};
    end

    idx = find(strcmp({meta.fecha}, date) & [meta.planta] == pot, 1);
    if isempty(idx)
        error('Pot %d, date %s, side %s not found.', pot, date, side);
    end

    row = meta(idx).fila;
    col  = meta(idx).columna;

    img = zeros(nRows, nCols, 3);

    if ~iscell(info_table.date)
        info_table.date = cellstr(info_table.date);
    end
    if ~iscell(info_table.leaf)
        info_table.leaf = cellstr(info_table.leaf);
    end

    for i = 1:3
        pix = pixin{i}{row, col};
        if isempty(pix), continue; end

        leaf_str = sprintf('Leaf%d%s', i, side);
        idxInfo = find(info_table.pot == pot & ...
                       strcmp(info_table.date, date) & ...
                       strcmp(info_table.leaf, leaf_str));
        if isempty(idxInfo), continue; end

        T2_vals = info_table.T2_all{idxInfo};

        maskGreen = false(nRows, nCols);
        maskRed  = false(nRows, nCols);
        maskGreen(pix) = T2_vals <= T2_limit;
        maskRed(pix)  = T2_vals >  T2_limit;

        img(:,:,2) = img(:,:,2) + maskGreen;   % green
        img(:,:,1) = img(:,:,1) + maskRed;    % red
    end
end
