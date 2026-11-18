function T = build_leaf_table(result, cfg)
%BUILD_LEAF_TABLE  One row per leaf x date observation
%
%   T = build_leaf_table(result, cfg)
%
%   Goes through 'result' (one entry per pot x date) and locates, inside
%   result{i}.data, the block of pixels of each leaf. Leaves are stored
%   in the order: image G (left side) leaves 1, 2, 3, then image D (right
%   side) leaves 1, 2, 3; leaves with 0 pixels are skipped.
%
%   Columns of T:
%     i_res      index in 'result'
%     pot     pot number
%     date      date 'yyyymmdd'
%     dpi        days post-inoculation
%     leaf       leaf label 'Leaf1G' ... 'Leaf3D'
%     side       image side, 'G' (left) or 'D' (right)
%     leaf_id    unique leaf identifier, e.g. 'P25_Leaf1G'
%     image_id  image identifier, e.g. 'P25_20210311_G'
%     genotype   'EPO67' or 'EPO68'
%     inoculated  true (inoculated) / false (mock-inoculated)
%     npix, row_start, row_end  number of pixels and rows within result{i}.data
%
%   Pots not assigned to any group in cfg are removed (with a warning), and
%   leaves listed in cfg.exclude are removed.

sides = {'G','G','G','D','D','D'};
nums  = [1 2 3 1 2 3];

rows = cell(0, 8);
for i = 1:numel(result)
    e = result{i};
    if ~ismember(e.date, cfg.dates)
        continue;                      % dates outside the analysis (30-31/03)
    end
    row_start = 1;
    for h = 1:6
        field_name = sprintf('pixelsLeaf%d%s', nums(h), sides{h});
        n = 0;
        if isfield(e, field_name), n = e.(field_name); end
        if n <= 0, continue; end
        row_end = row_start + n - 1;
        rows(end+1, :) = {i, e.pot, e.date, ...
            sprintf('Leaf%d%s', nums(h), sides{h}), sides{h}, n, row_start, row_end}; %#ok<AGROW>
        row_start = row_end + 1;
    end
    % Check: the leaf blocks must cover 'data' exactly
    if row_start - 1 ~= size(e.data, 1)
        warning('Pot %d, date %s: %d pixels in leaves, %d in data.', ...
            e.pot, e.date, row_start - 1, size(e.data, 1));
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'i_res','pot','date','leaf','side','npix','row_start','row_end'});

% --- dpi ---
[~, loc] = ismember(T.date, cfg.dates);
T.dpi = cfg.dpi(loc)';

% --- leaf and image identifiers ---
T.leaf_id   = strcat('P', string(T.pot), '_', string(T.leaf));
T.image_id = strcat('P', string(T.pot), '_', string(T.date), '_', string(T.side));
T.leaf_id   = cellstr(T.leaf_id);
T.image_id = cellstr(T.image_id);

% --- genotype and treatment ---
n = height(T);
T.genotype  = repmat({''}, n, 1);
T.inoculated = false(n, 1);
for r = 1:n
    p = T.pot(r);
    if ismember(p, [cfg.inoc_EPO67 cfg.mock_EPO67])
        T.genotype{r} = 'EPO67';
    elseif ismember(p, [cfg.inoc_EPO68 cfg.mock_EPO68])
        T.genotype{r} = 'EPO68';
    end
    T.inoculated(r) = ismember(p, [cfg.inoc_EPO67 cfg.inoc_EPO68]);
end
no_group = cellfun(@isempty, T.genotype);
if any(no_group)
    fprintf('Pots not used in this run (excluded from the reference): %s\n', ...
        mat2str(unique(T.pot(no_group))'));
    T(no_group, :) = [];
end

% --- excluded leaves ---
for k = 1:size(cfg.exclude, 1)
    to_remove = T.pot == cfg.exclude{k, 1} & strcmp(T.leaf, cfg.exclude{k, 2});
    T(to_remove, :) = [];
end
end
