function run_preprocessing(varargin)
%RUN_PREPROCESSING  From the raw VNIR images of the database to alldataandmeta.mat
%
%   run_preprocessing                 all images of the 20 pots (edit the default
%                                     folders below, or pass them as options)
%   run_preprocessing(Name, Value)    options:
%     'raw_dir'      downloaded images: raw_dir/<yyyymmdd>/<pot>_<side>/...
%     'work_dir'     results per image (small files)
%     'pots'         pots processed (default: the 20 pots; [] = all pots)
%     'dates'        dates processed, e.g. {'20210316'} (default: all, except April)
%     'output_file'  default: alldataandmeta.mat in the code folder
%
% One function for the whole preprocessing. For each image of the selected pots
% and dates it:
%   1. reads the raw HySpex VNIR image (1024 x 1370 x 216, uint16; one image,
%      20210312 29_G, has 1349 lines: missing columns are left empty) with hsread;
%   2. converts it to reflectance with the Spectralon reference (mean of the
%      first 50 columns, line by line);
%   3. crops the inoculated area (rect = [0 110 1370 780], 781 x 1370 pixels);
%   4. reads the leaf masks of the database (<pot>_<side>_leaf1/2/3.csv: x =
%      column, y = row in the cropped image), keeping the leaf numbers;
%   5. builds the cropped cube with the leaf spectra and zeros elsewhere,
%      resizes it with imresize (factor 0.4, bicubic) to 313 x 548 x 216 and
%      maps the leaf pixels to the resized image (x 0.4, rounded);
%   6. stores the spectra and pixel indices of each leaf.
% At the end, all images are organised into alldataandmeta.mat (code folder),
% the input of run_analysis.m.
%
% The intermediate full-size cubes are not written to disk. The results of each
% image are saved in work_dir/per_image/ (a few MB each), so the script can be
% stopped and restarted: images already processed are not processed again.
%
% Requirements: Image Processing Toolbox; hsread (preprocessing/hsi_toolbox/).
% Runs on MATLAB R2015a or later.

%% Settings (default folders: edit, or pass them as options)
ip = inputParser;
ip.addParameter('raw_dir',  'C:\data\STB\raw_VNIR');
ip.addParameter('work_dir', 'C:\data\STB\preprocessing');
% the 20 pots imaged from -1 to 17 dpi (the analysis then excludes the
% mock-inoculated pots flagged by the pre-inoculation check)
ip.addParameter('pots', [25 26 27 28 29 32 33 34 35 36 37 38 41 42 43 44 45 46 47 48]);
ip.addParameter('dates', {});
ip.addParameter('output_file', '');
ip.parse(varargin{:});
raw_dir = ip.Results.raw_dir; work_dir = ip.Results.work_dir; pots_used = ip.Results.pots;

%% Constants (as in the original scripts)
ncol = 50;                   % Spectralon columns
rect = [0 110 1370 780];     % crop [left top width height] -> 781 x 1370
Y = 781; X = 1370; Z = 216;  % cropped cube size
f = 0.4;                     % resizing factor

here     = fileparts(mfilename('fullpath'));
code_dir = fileparts(here);
addpath(fullfile(here, 'hsi_toolbox'));                      % hsread
img_dir = fullfile(work_dir, 'per_image');
if ~exist(img_dir, 'dir'), mkdir(img_dir); end

%% Dates: folders named yyyymmdd; April (after 17 dpi) is not used
d = dir(fullfile(raw_dir, '20*'));
d = d([d.isdir]);
dates = sort({d.name});
dates = dates(~cellfun(@(n) strcmp(n(5:6), '04'), dates));
if ~isempty(ip.Results.dates), dates = intersect(dates, cellstr(ip.Results.dates)); end

%% 1-6. Process each image
for i = 1:numel(dates)
    files = rdir(fullfile(raw_dir, dates{i}), '*.hyspex');   % also in subfolders (<pot>_<side>/)
    done_here = {};
    for j = 1:numel(files)
        [~, name] = fileparts(files(j).name);                      % e.g. 33_D
        pot = str2double(strtok(name, '_'));
        if ~isempty(pots_used) && ~ismember(pot, pots_used), continue; end
        if ismember(name, done_here)
            warning('%s %s found twice in the date folder; second copy ignored (%s)', dates{i}, name, files(j).folder);
            continue
        end
        csvbase = fullfile(files(j).folder, name);
        if isempty(dir([csvbase '_leaf*.csv']))
            fprintf('%s %s: no leaf mask in the database, image skipped\n', dates{i}, name);
            continue
        end
        done_here{end+1} = name; %#ok<SAGROW>
        out_file = fullfile(img_dir, sprintf('img_%s_%s_resized.mat', dates{i}, name(1:4)));
        if exist(out_file, 'file'), continue; end                  % already processed (restart)
        fprintf('%s %s\n', dates{i}, name);

        % 1-2. Raw image and reflectance
        hs = double(hsread(fullfile(files(j).folder, files(j).name))) / 65536;
        if size(hs, 1) ~= 1024 || size(hs, 3) ~= Z
            warning('%s %s: unexpected image size %s, image skipped', dates{i}, name, mat2str(size(hs)));
            continue
        end
        hs_R = bsxfun(@rdivide, hs, mean(hs(:, 1:ncol, :), 2));
        clear hs

        % 3-4. Leaf masks of the database in the cropped image
        sz_crop = size(imcrop(hs_R(:,:,1), rect));
        L = zeros(sz_crop);
        for k = 1:3
            fcsv = sprintf('%s_leaf%d.csv', csvbase, k);
            if exist(fcsv, 'file') == 2
                P = csvread(fcsv);                                 % x (column), y (row)
                out = P(:,1) < 1 | P(:,1) > sz_crop(2) | P(:,2) < 1 | P(:,2) > sz_crop(1);
                if any(out)
                    warning('%s %s leaf %d: %d mask pixels outside the image (%d x %d) ignored', ...
                        dates{i}, name, k, sum(out), sz_crop(1), sz_crop(2));
                    P = P(~out, :);
                end
                L(sub2ind(sz_crop, P(:,2), P(:,1))) = k;
            end
        end
        L = L(1:min(Y, end), 1:min(X, end));                       % cube bounds, as in the original

        % 5. Cropped cube with the leaf spectra (zeros elsewhere), resized
        M = L > 0;
        [ny, nx] = size(L);
        cube = zeros(Y, X, Z);
        for z = 1:Z
            band = imcrop(hs_R(:,:,z), rect);
            band = band(1:ny, 1:nx);
            t = zeros(ny, nx);
            t(M) = band(M);                                        % leaf pixels only
            cube(1:ny, 1:nx, z) = t;
        end
        clear hs_R
        cube = imresize(cube, f);
        sz = size(cube);
        Xr = reshape(cube, [], Z);
        clear cube

        % leaf pixel coordinates (x = column, y = row) and their indices in the resized image
        getResizedIdx = @(pix) unique(sub2ind(sz(1:2), ...
            min(max(round(pix(:,2) * f), 1), sz(1)), ...
            min(max(round(pix(:,1) * f), 1), sz(2))));
        img_mask = zeros(sz(1:2));
        pix = cell(1, 3); idx = cell(1, 3); Zs = cell(1, 3);
        for k = 1:3
            [rr, cc] = find(L == k);
            pix{k} = [cc rr];
            for m = 1:size(pix{k}, 1)
                yy = round(pix{k}(m,2) * f); xx = round(pix{k}(m,1) * f);
                if xx >= 1 && xx <= sz(2) && yy >= 1 && yy <= sz(1), img_mask(yy, xx) = 1; end
            end
            if isempty(pix{k}), idx{k} = []; else, idx{k} = getResizedIdx(pix{k}); end
            Zs{k} = Xr(idx{k}, :);
        end
        pixin_all = find(img_mask(:) > 0);

        % 6. Results of this image
        Z1 = Zs{1}; Z2 = Zs{2}; Z3 = Zs{3};
        pixin1_resized = idx{1}; pixin2_resized = idx{2}; pixin3_resized = idx{3};
        pixout_all = setdiff(1:numel(img_mask), pixin_all)';
        save(out_file, 'Z1', 'Z2', 'Z3', 'pixin1_resized', 'pixin2_resized', 'pixin3_resized', 'pixout_all');
    end
end

%% Organise all images into alldataandmeta.mat
% Same variables, order and metadata fields as the deposited file:
% cell arrays {date, image}, with the images of each side sorted by file name.
% Field names of the metadata (kept from the original file): fecha (date),
% tipo (side), fila (date index), columna (image index), planta (pot),
% filename, variables_cargadas (variables found), faltantes (variables missing).
expected = {'Z1','Z2','Z3','pixin1_resized','pixin2_resized','pixin3_resized','pixout_all'};
D_Z1 = {}; D_Z2 = {}; D_Z3 = {}; D_pixin1 = {}; D_pixin2 = {}; D_pixin3 = {}; D_pixout = {};
G_Z1 = {}; G_Z2 = {}; G_Z3 = {}; G_pixin1 = {}; G_pixin2 = {}; G_pixin3 = {}; G_pixout = {};
metadataD = struct('fecha', {}, 'tipo', {}, 'fila', {}, 'planta', {}, ...
                   'filename', {}, 'variables_cargadas', {}, 'faltantes', {});
metadataG = metadataD;

has_img = false(size(dates));
for i = 1:numel(dates)
    fl = dir(fullfile(img_dir, ['img_' dates{i} '_*_resized.mat']));
    pf = cellfun(@(n) str2double(regexp(n, '^img_\d{8}_(\d+)_', 'tokens', 'once')), {fl.name});
    has_img(i) = any(isempty(pots_used) | ismember(pf, pots_used));
end
dateFolderNames = dates(has_img);
for fi = 1:numel(dateFolderNames)
    for side = ['D', 'G']
        files = dir(fullfile(img_dir, ['img_' dateFolderNames{fi} '_*_' side '_resized.mat']));
        [~, o] = sort({files.name}); files = files(o);
        if ~isempty(pots_used)    % only the selected pots (work_dir may hold others)
            pf = cellfun(@(n) str2double(regexp(n, ['(\d+)_' side '_resized'], 'tokens', 'once')), {files.name});
            files = files(ismember(pf, pots_used));
        end
        for n = 1:numel(files)
            data = load(fullfile(img_dir, files(n).name));
            varNames = fieldnames(data);
            tokens = regexp(files(n).name, ['(\d+)_' side], 'tokens');
            entry = struct();
            entry.fecha = dateFolderNames{fi};
            entry.tipo = side;
            entry.fila = fi;
            entry.columna = n;
            entry.planta = str2double(tokens{1}{1});
            entry.filename = files(n).name;
            entry.variables_cargadas = intersect(expected, varNames);
            entry.faltantes = expected(~ismember(expected, varNames));
            if side == 'D'
                D_Z1{fi,n} = data.Z1; D_Z2{fi,n} = data.Z2; D_Z3{fi,n} = data.Z3;
                D_pixin1{fi,n} = data.pixin1_resized; D_pixin2{fi,n} = data.pixin2_resized;
                D_pixin3{fi,n} = data.pixin3_resized; D_pixout{fi,n} = data.pixout_all;
                if isempty(metadataD), metadataD = entry; else, metadataD(end+1) = entry; end %#ok<SAGROW>
            else
                G_Z1{fi,n} = data.Z1; G_Z2{fi,n} = data.Z2; G_Z3{fi,n} = data.Z3;
                G_pixin1{fi,n} = data.pixin1_resized; G_pixin2{fi,n} = data.pixin2_resized;
                G_pixin3{fi,n} = data.pixin3_resized; G_pixout{fi,n} = data.pixout_all;
                if isempty(metadataG), metadataG = entry; else, metadataG(end+1) = entry; end %#ok<SAGROW>
            end
        end
    end
end
output_file = ip.Results.output_file;
if isempty(output_file), output_file = fullfile(code_dir, 'alldataandmeta.mat'); end
save(output_file, ...
    'D_Z1','D_Z2','D_Z3','D_pixin1','D_pixin2','D_pixin3','D_pixout', ...
    'G_Z1','G_Z2','G_Z3','G_pixin1','G_pixin2','G_pixin3','G_pixout', ...
    'dateFolderNames','metadataD','metadataG', '-v7.3');
fprintf('Saved %s: %d dates, %d images (D) + %d images (G)\n', ...
    output_file, numel(dateFolderNames), numel(metadataD), numel(metadataG));

%% Wavelengths of the 216 bands (for the plots), added to soptfullmodel.mat if missing
sopt_file = fullfile(code_dir, 'soptfullmodel.mat');
if exist(sopt_file, 'file') && ~ismember('wavelengths', who('-file', sopt_file))
    h = rdir(raw_dir, '*.hdr');
    txt = fileread(fullfile(h(1).folder, h(1).name));
    tok = regexp(txt, 'wavelength\s*=\s*\{([^}]*)\}', 'tokens', 'once', 'ignorecase');
    wavelengths = str2double(strsplit(strtrim(tok{1}), ','));
    wavelengths = wavelengths(~isnan(wavelengths));
    save(sopt_file, 'wavelengths', '-append');
    fprintf('wavelengths added to %s\n', sopt_file);
end
end

function out = rdir(root, pattern)
% Files matching pattern in root and in all its subfolders (fields name and
% folder). Written for MATLAB versions without the recursive dir ('**',
% R2016b) and without the field 'folder' of dir.
out = struct('name', {}, 'folder', {});
d = dir(fullfile(root, pattern));
for i = 1:numel(d)
    if ~d(i).isdir, out(end+1) = struct('name', d(i).name, 'folder', root); end %#ok<AGROW>
end
s = dir(root);
for i = 1:numel(s)
    if s(i).isdir && ~any(strcmp(s(i).name, {'.', '..'}))
        out = [out, rdir(fullfile(root, s(i).name), pattern)]; %#ok<AGROW>
    end
end
end
