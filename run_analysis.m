function run_analysis(varargin)
%RUN_ANALYSIS  Complete analysis of the article (MCR-ALS + MSPC, leave-one-pot-out)
%
%   run_analysis                  runs all the steps below
%   run_analysis('steps', {'figures'})   runs only some steps (they read the
%                                 results saved by the previous ones)
%   run_analysis(Name, Value)     changes an option (see "Options")
%
% Input file (code folder): alldataandmeta.mat (preprocessing/run_preprocessing.m).
% The profiles of the full multiset MCR-ALS model and the wavelengths are
% written in functions/full_model_profiles.m.
% All results are written to the output folder (results/ by default).
%
% Steps, in this order:
%   'data'         leaf spectra per pot x date (revision_data.mat)
%   'check'        pre-inoculation check of the mock-inoculated pots, applied
%                  iteratively: the pot with the most significant excess of
%                  alarms at -1 and 0 dpi is removed and the check repeated
%                  until no pot is flagged. The first iteration (all pots) is
%                  also the sensitivity analysis with all 20 pots.
%   'main'         design (Table 1), number of components and stability
%                  (Figs. S1-S3), leave-one-pot-out validation (Table 3,
%                  Table S2), mixed models (Table 2, Table S5), pixel
%                  assignment (Section 3.3), references updated over time
%                  (Table S4) and metadata files (metadata/)
%   'sensitivity'  three PCs with T2 and SPE, SIMPLISMA initial estimates in
%                  each fold, non-negative projection; summary (Table S3)
%   'figures'      Figures 2-6
%
% Options (Name, Value):
%   'steps'            cell array of steps (default: all)
%   'output_folder'    default 'results'
%   'n_boot'           bootstrap resamples of pots (default 2000)
%   'n_rep_stability'  random initialisations for Fig. S3 (default 20)
%   'max_pix_mcr'      pixels per leaf and date used to fit MCR-ALS (default 200)
%   'check_reference'  true: run the iterative pre-inoculation check;
%                      false: use the exclusions of the article (pots 36, 46)
%   'seed'             random seed (default 2026)
%
% Helper functions are in functions/. Runs on MATLAB R2020a or later with the
% Statistics and Machine Learning Toolbox.

%% Options
ip = inputParser;
ip.addParameter('steps', {'data','check','main','sensitivity','figures'});
ip.addParameter('output_folder', 'results');
ip.addParameter('n_boot', 2000);
ip.addParameter('n_rep_stability', 20);
ip.addParameter('max_pix_mcr', 200);
ip.addParameter('check_reference', true);
ip.addParameter('seed', 2026);
ip.parse(varargin{:});
opt = ip.Results;
steps = cellstr(opt.steps);

cfg = default_config(opt);
remove_dataset_conflicts();
t_start = tic;

%% Data
if ismember('data', steps)
    step_organise_data(cfg);
end

%% Pre-inoculation check of the reference (iterative)
excl_file = fullfile(cfg.output_folder, 'excluded_pots.mat');
if ismember('check', steps) && opt.check_reference
    excluded = step_reference_check_iterative(cfg);
    save(excl_file, 'excluded');
elseif exist(excl_file, 'file') == 2
    load(excl_file, 'excluded');
else
    excluded = struct('pot', {36, 46}, 'genotype', {'EPO67', 'EPO67'}, 'iteration', {1, 2});
end
cfg = apply_exclusions(cfg, excluded);
fprintf('\nMain analysis: excluded mock-inoculated pots %s\n', mat2str([excluded.pot]));

%% Main analysis
if ismember('main', steps)
    step_design_table(cfg);
    step_mcr_selection_stability(cfg, opt.n_rep_stability);
    step_leave_one_pot_out(cfg);
    step_reference_check(cfg);          % final check (no pot should be flagged)
    step_mixed_models(cfg);
    step_pixels(cfg);
    step_time_updated_reference(cfg);
    step_export_metadata(cfg, excluded);
end

%% Sensitivity analyses (one option changed at a time)
if ismember('sensitivity', steps)
    c = cfg; c.mspc.k = [];               step_leave_one_pot_out(c);   % three PCs (>= 95 %), T2 and SPE
    c = cfg; c.mcr_init = 'simplisma';    step_leave_one_pot_out(c);   % SIMPLISMA in each fold
    c = cfg; c.projection = 'nnls';       step_leave_one_pot_out(c);   % non-negative projection
    step_sensitivity_summary(cfg);
end

%% Figures
if ismember('figures', steps)
    fig2_trajectories(cfg);
    fig3_loadings(cfg);
    fig4_T2(cfg);
    fig5_pixels(cfg);
    fig6_maps(cfg);
end
fprintf('\nDone in %.1f h. Results in %s\n', toc(t_start) / 3600, cfg.output_folder);
end


%% ========================================================================
%  CONFIGURATION
%  ========================================================================
function cfg = default_config(opt)
cfg.code_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(cfg.code_dir, 'functions'));
cfg.alldata_file = fullfile(cfg.code_dir, 'alldataandmeta.mat');
cfg.data_file    = fullfile(cfg.code_dir, 'revision_data.mat');   % written by step 'data'
if startsWith(opt.output_folder, {'/', '\'}) || contains(opt.output_folder, ':')
    cfg.output_folder = opt.output_folder;                         % absolute path
else
    cfg.output_folder = fullfile(cfg.code_dir, opt.output_folder); % relative to the code folder
end
if ~exist(cfg.output_folder, 'dir'), mkdir(cfg.output_folder); end

% Wavelengths: bands 25-130 of the 216 VNIR bands (476-760 nm)
cfg.idx_lambda = 25:130;
[~, w] = full_model_profiles();
cfg.lambda_nm = w(cfg.idx_lambda);

% Dates and dpi (inoculation on 12/03/2021 = 0 dpi)
cfg.dates = {'20210311','20210312','20210316','20210317','20210318', ...
             '20210319','20210320','20210322','20210323','20210324', ...
             '20210325','20210326','20210327','20210328','20210329'};
cfg.dpi   = [-1 0 4 5 6 7 8 10 11 12 13 14 15 16 17];

% Design: the 20 pots imaged from -1 to 17 dpi (numbers of the database,
% doi:10.57745/WVP0FJ). Mock-inoculated pots are removed by the
% pre-inoculation check (step 'check').
cfg.inoc_EPO67 = [32 33 43 44 45];
cfg.inoc_EPO68 = [25 26 27 37 38];
cfg.mock_EPO67 = [34 35 36 46 47 48];
cfg.mock_EPO68 = [28 29 41 42];
cfg.exclude = {};                 % excluded leaves, e.g. {46, 'Leaf2G'}
% Visual confidence level of the infection status of each leaf (database file
% Data_description.csv: 1 = high, 2 = low). Leaves with a low confidence level
% are not used, because their status (healthy or infected) is uncertain: a
% doubtful mock-inoculated leaf would contaminate the healthy reference and a
% doubtful inoculated leaf would bias the detection metrics.
% The file is downloaded from the database and placed in the code folder
% (or in database/).
cfg.confidence_file = fullfile(cfg.code_dir, 'Data_description.csv');
if exist(cfg.confidence_file, 'file') ~= 2
    cfg.confidence_file = fullfile(cfg.code_dir, 'database', 'Data_description.csv');
end
if exist(cfg.confidence_file, 'file') ~= 2
    error(['Data_description.csv not found. Download it from the database ' ...
           '(https://doi.org/10.57745/WVP0FJ) and place it in the code folder.']);
end

% MCR-ALS
cfg.nc               = 4;
cfg.max_pix_mcr      = opt.max_pix_mcr;   % pixels per leaf and date used to fit MCR-ALS
cfg.n_simplisma      = 20000;             % candidate pixels for SIMPLISMA
cfg.simplisma_offset = 10;                % SIMPLISMA noise level (% of the maximum mean)
cfg.als.max_iter     = 200;
cfg.als.tol          = 0.1;               % relative change in lack of fit < 0.1 %

% Validation (main analysis; the sensitivity analyses change one option)
cfg.mcr_mock_only = false;        % true: MCR-ALS on mock-inoculated calibration leaves only
cfg.mcr_init   = 'full_model';    % 'full_model' or 'simplisma' (initial estimates in each fold)
cfg.projection = 'pinv';          % 'pinv' (least squares) or 'nnls' (non-negative)
cfg.seed   = opt.seed;
cfg.n_boot = opt.n_boot;

% MSPC
cfg.mspc.alpha   = 0.05;
cfg.mspc.var_pca = 0.95;          % used when k is empty
cfg.mspc.k       = 4;             % number of PCs ([] = minimum explaining >= 95 %)
cfg.mspc.limit   = 'F';           % phase-II F limit for new observations
cfg.mspc.use_Q   = true;          % alarm if T2 or SPE above its limit (SPE = 0 with 4 PCs)
cfg.mspc.ref_by_genotype = false;

% Pixel maps
cfg.pix.percentile = 99;          % pixel T2 limit: 99th percentile of reference pixels

% Quality check of the leaf masks: a leaf x date observation is excluded when
% fewer than 60 % of its pixels reflect more than 0.15 at 760 nm, i.e. when
% its mask lies mostly on the black background instead of on the leaf (leaves
% reflect about 0.3-0.5 at 760 nm, still about 0.2 when senescent).
cfg.qc.band            = 130;     % band 130 of the 216 VNIR bands = 760 nm
cfg.qc.min_reflectance = 0.15;
cfg.qc.min_fraction    = 0.60;

% Figure font
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');
set(groot, 'defaultColorbarFontName', 'Times New Roman');
end

function cfg = apply_exclusions(cfg, excluded)
for e = excluded(:)'
    cfg.mock_EPO67 = setdiff(cfg.mock_EPO67, e.pot, 'stable');
    cfg.mock_EPO68 = setdiff(cfg.mock_EPO68, e.pot, 'stable');
end
end

function C = read_confidence(file)
% Pot, side, leaf and confidence level from the database file
% Data_description.csv (semicolon-separated, notes in extra columns).
txt = fileread(file);
tok = regexp(txt, '^(\d+);([GD]);(\d);EPO_\d+;[^;]*;([12])', 'tokens', 'lineanchors');
tok = vertcat(tok{:});
C = table(str2double(tok(:, 1)), tok(:, 2), str2double(tok(:, 3)), str2double(tok(:, 4)), ...
    'VariableNames', {'pot','side','leaf','confidence'});
end

function remove_dataset_conflicts()
% Some toolboxes (e.g. PLS_Toolbox) define their own 'dataset' class, which
% shadows MATLAB's and makes fitlme fail. Their folders are removed from the
% path for this session only.
stats_dataset = {fullfile(matlabroot, 'toolbox', 'stats'), fullfile(matlabroot, 'toolbox', 'shared'), ...
                 fullfile(matlabroot, 'toolbox', 'matlab')};
cand = [which('dataset', '-all'); which('table2dataset', '-all')];
for i = 1:numel(cand)
    p = cand{i};
    if ~any(startsWith(p, stats_dataset)) && ~contains(p, 'built-in')
        folder = fileparts(p);
        while contains(folder, '@'), folder = fileparts(folder); end
        path_list = strsplit(path, pathsep);
        to_remove = path_list(startsWith(path_list, folder, 'IgnoreCase', true));
        if ~isempty(to_remove)
            warning('Removing from the path (conflicts with fitlme): %s (%d folders). If fitlme still fails, restart MATLAB.', ...
                folder, numel(to_remove));
            rmpath(strjoin(to_remove, pathsep));
        end
    end
end
end


%% ========================================================================
%  STEP 'data': leaf spectra per pot x date
%  ========================================================================
function step_organise_data(cfg)
% One entry per pot x date with the pixel spectra of its leaves (side G
% leaves 1-3, then side D leaves 1-3) and the number of pixels of each leaf.
% The full-model profiles (soptTots) are those of full_model_profiles.m.
fprintf('\n===== Organising the data =====\n');
rng(cfg.seed);
load(cfg.alldata_file, 'metadataG', 'metadataD', ...
     'G_Z1', 'G_Z2', 'G_Z3', 'D_Z1', 'D_Z2', 'D_Z3', ...
     'G_pixin1', 'G_pixin2', 'G_pixin3', 'D_pixin1', 'D_pixin2', 'D_pixin3', ...
     'dateFolderNames');
pots  = [28 29 34 35 36 41 42 46 47 48 25 26 27 32 33 37 38 43 44 45];
dates = dateFolderNames;

result = {};
n = 0;
for pot = pots
    for f = 1:numel(dates)
        [dataG, leavesG] = get_leaf_data('G', dates{f}, pot, metadataG, G_Z1, G_Z2, G_Z3);
        [dataD, leavesD] = get_leaf_data('D', dates{f}, pot, metadataD, D_Z1, D_Z2, D_Z3);
        if isempty(dataG) && isempty(dataD), continue; end
        n = n + 1;
        result{n} = struct('pot', pot, 'date', dates{f}, ...
            'data', [dataG; dataD], 'leavesG', leavesG, 'leavesD', leavesD); %#ok<AGROW>
    end
end

pixG = {G_pixin1, G_pixin2, G_pixin3};
pixD = {D_pixin1, D_pixin2, D_pixin3};
npix_all = 0;
for i = 1:numel(result)
    e = result{i};
    for leaf = 1:3
        result{i}.(sprintf('pixelsLeaf%dG', leaf)) = 0;
        result{i}.(sprintf('pixelsLeaf%dD', leaf)) = 0;
    end
    iG = find(strcmp({metadataG.fecha}, e.date) & [metadataG.planta] == e.pot, 1);
    iD = find(strcmp({metadataD.fecha}, e.date) & [metadataD.planta] == e.pot, 1);
    for leaf = 1:3
        if ~isempty(iG) && ismember(leaf, e.leavesG)
            result{i}.(sprintf('pixelsLeaf%dG', leaf)) = ...
                numel(pixG{leaf}{metadataG(iG).fila, metadataG(iG).columna});
        end
        if ~isempty(iD) && ismember(leaf, e.leavesD)
            result{i}.(sprintf('pixelsLeaf%dD', leaf)) = ...
                numel(pixD{leaf}{metadataD(iD).fila, metadataD(iD).columna});
        end
    end
    npix = sum(cellfun(@(c) result{i}.(c), ...
        {'pixelsLeaf1G','pixelsLeaf2G','pixelsLeaf3G','pixelsLeaf1D','pixelsLeaf2D','pixelsLeaf3D'}));
    if npix ~= size(e.data, 1)
        warning('Pot %d, date %s: %d pixels in masks, %d spectra.', e.pot, e.date, npix, size(e.data, 1));
    end
    npix_all = npix_all + size(e.data, 1);
end
fprintf('%d pot x date entries, %d pixel spectra.\n', numel(result), npix_all);

soptTots = full_model_profiles();
save(cfg.data_file, 'result', 'soptTots', '-v7.3');
fprintf('Saved %s\n', cfg.data_file);
end


%% ========================================================================
%  STEP 'check': pre-inoculation check of the reference, iterative
%  ========================================================================
function excluded = step_reference_check_iterative(cfg)
% A mock-inoculated pot is flagged when the proportion of its observations
% at -1 and 0 dpi above the control limit is significantly higher than
% alpha (one-sided binomial test, p < 0.05). An anomalous pot widens the
% reference and can mask another one, so the most significant pot is
% removed and the check is repeated until no pot is flagged. Uses only
% data acquired before inoculation.
excluded = struct('pot', {}, 'genotype', {}, 'iteration', {});
base = cfg.output_folder;
it = 0;
while true
    it = it + 1;
    c = apply_exclusions(cfg, excluded);
    if it == 1
        c.output_folder = fullfile(base, 'all_20_pots');   % also a sensitivity analysis
    else
        c.output_folder = fullfile(base, 'reference_check', sprintf('iteration_%d', it));
    end
    if ~exist(c.output_folder, 'dir'), mkdir(c.output_folder); end
    fprintf('\n===== Pre-inoculation check, iteration %d (excluded: %s) =====\n', ...
        it, mat2str([excluded.pot]));
    step_design_table(c);
    step_leave_one_pot_out(c);
    Chk = step_reference_check(c);
    if it == 1
        step_mixed_models(c);           % divergence dates with all 20 pots
    end
    if ~any(Chk.exclude), break; end
    k = find(Chk.exclude, 1);            % Chk is sorted by p-value
    excluded(end+1) = struct('pot', Chk.pot_id(k), 'genotype', Chk.genotype{k}, 'iteration', it); %#ok<AGROW>
    fprintf('Pot %d removed from the reference.\n', Chk.pot_id(k));
end
end

function Chk = step_reference_check(cfg)
run_mode = mode_name(cfg);
load(fullfile(cfg.output_folder, sprintf('validation_%s.mat', run_mode)), 'R');
alpha = cfg.mspc.alpha;
pre = R(~R.inoculated & R.dpi <= 0, :);
[g, pot_id, genotype] = findgroups(pre.pot, pre.genotype);
n_pre      = splitapply(@numel, pre.alarm, g);
alarms     = splitapply(@sum, pre.alarm, g);
p_binomial = 1 - binocdf(alarms - 1, n_pre, alpha);   % P(X >= alarms)
exclude    = p_binomial < 0.05;
Chk = table(pot_id, genotype, n_pre, alarms, round(100 * alarms ./ n_pre, 1), ...
    p_binomial, exclude, 'VariableNames', ...
    {'pot_id','genotype','n_obs_pre','alarms_pre','pct_alarms_pre','p_binomial','exclude'});
Chk = sortrows(Chk, 'p_binomial');
disp(Chk);
if any(Chk.exclude)
    fprintf('Flagged pot(s): %s\n', mat2str(Chk.pot_id(Chk.exclude)'));
else
    fprintf('No mock-inoculated pot differs from the reference before inoculation.\n');
end
writetable(Chk, fullfile(cfg.output_folder, sprintf('Reference_check_%s.xlsx', run_mode)));
end


%% ========================================================================
%  STEP 'main'
%  ========================================================================
function step_design_table(cfg)
% Leaf x date table (used by all the other steps) and Table 1.
load(cfg.data_file, 'result');
T = build_leaf_table(result, cfg);
T.treatment = repmat({'mock-inoculated'}, height(T), 1);
T.treatment(T.inoculated) = {'inoculated'};

% Leaves with a low confidence level of their infection status
C = read_confidence(cfg.confidence_file);
key_T = strcat(string(T.pot), '_', string(T.leaf));
key_C = strcat(string(C.pot), '_Leaf', string(C.leaf), C.side);
low = ismember(key_T, key_C(C.confidence == 2));
LC = unique(T(low, {'pot','side','leaf','genotype','treatment'}));
LC.n_observations = zeros(height(LC), 1);
for r = 1:height(LC)
    LC.n_observations(r) = sum(low & T.pot == LC.pot(r) & strcmp(T.leaf, LC.leaf{r}));
end
writetable(LC, fullfile(cfg.output_folder, 'excluded_leaves_low_confidence.csv'));
fprintf('\n===== Confidence level of the infection status =====\n');
fprintf('%d leaves (%d observations) with a low confidence level excluded:\n', height(LC), sum(low));
disp(LC);
T(low, :) = [];

% Quality check of the leaf masks (cfg.qc): observations whose mask lies
% mostly on the background are excluded
frac = zeros(height(T), 1);
for r = 1:height(T)
    frac(r) = mean(result{T.i_res(r)}.data(T.row_start(r):T.row_end(r), cfg.qc.band) > cfg.qc.min_reflectance);
end
bad = frac < cfg.qc.min_fraction;
QC = T(bad, {'pot','side','leaf','date','dpi','genotype','treatment','npix'});
QC.pct_pixels_on_leaf = round(100 * frac(bad), 1);
writetable(QC, fullfile(cfg.output_folder, 'excluded_leaf_observations.csv'));
fprintf('\n===== Quality check of the leaf masks =====\n');
fprintf('%d of %d leaf x date observations excluded (fewer than %.0f %% of the pixels above %.2f at %.0f nm):\n', ...
    sum(bad), height(T), 100 * cfg.qc.min_fraction, cfg.qc.min_reflectance, cfg.lambda_nm(cfg.qc.band - cfg.idx_lambda(1) + 1));
disp(QC);
[fs, o] = sort(frac(~bad)); Tk = T(~bad, :);
fprintf('Lowest values among the observations kept: %s %%\n', mat2str(round(100 * fs(1:min(5, end))', 1)));
fprintf('  (%s)\n', strjoin(strcat(Tk.date(o(1:min(5, end))), {' pot '}, cellstr(num2str(Tk.pot(o(1:min(5, end))))), {' '}, Tk.leaf(o(1:min(5, end))))', '; '));
T(bad, :) = [];
save(fullfile(cfg.output_folder, 'leaf_table.mat'), 'T');

fprintf('\n===== Design =====\n');
fprintf('Pots %d, leaves %d, leaf x date observations %d, images %d (inoculated %d + mock %d), dates %d\n', ...
    numel(unique(T.pot)), numel(unique(T.leaf_id)), height(T), numel(unique(T.image_id)), ...
    numel(unique(T.image_id(T.inoculated))), numel(unique(T.image_id(~T.inoculated))), numel(unique(T.date)));
[~, ~, ih] = unique(T.leaf_id);
dpl = accumarray(ih, 1);
fprintf('Dates per leaf: min %d, median %d, max %d\n', min(dpl), median(dpl), max(dpl));

[g, geno, treat, dpi] = findgroups(T.genotype, T.treatment, T.dpi);
Table1 = table(geno, treat, dpi, splitapply(@(x) numel(unique(x)), T.pot, g), ...
    splitapply(@numel, T.leaf_id, g), splitapply(@(x) numel(unique(x)), T.image_id, g), ...
    'VariableNames', {'Genotype','Treatment','dpi','n_pots','n_leaves','n_images'});
Table1 = sortrows(Table1, {'Genotype','Treatment','dpi'});
tmp = Table1(:, {'Genotype','Treatment','dpi','n_leaves'});
tmp.dpi = categorical(tmp.dpi, cfg.dpi, strrep(compose('dpi_%d', cfg.dpi), '-', 'm'));
Table1_wide = unstack(tmp, 'n_leaves', 'dpi');
[g2, geno2, treat2] = findgroups(T.genotype, T.treatment);
Groups = table(geno2, treat2, splitapply(@(x) {mat2str(unique(x)')}, T.pot, g2), ...
    splitapply(@(x) numel(unique(x)), T.leaf_id, g2), ...
    'VariableNames', {'Genotype','Treatment','Pots','n_leaves'});
disp(Groups);
fx = fullfile(cfg.output_folder, 'Table1_design.xlsx');
writetable(Table1,      fx, 'Sheet', 'long_format');
writetable(Table1_wide, fx, 'Sheet', 'leaves_per_dpi');
writetable(Groups,      fx, 'Sheet', 'groups');
end

function step_mcr_selection_stability(cfg, n_rep)
% Singular values (Fig. S1), MCR-ALS with 3-5 components (Fig. S2) and
% stability against initialisation (Fig. S3). Descriptive only (no
% detection), so all the data are used.
fprintf('\n===== MCR-ALS: number of components and stability =====\n');
load(cfg.data_file, 'result', 'soptTots');
load(fullfile(cfg.output_folder, 'leaf_table.mat'), 'T');
rng(cfg.seed);
lambda_nm = cfg.lambda_nm;

[D, ~] = extract_pixels(result, T, 1:height(T), cfg.idx_lambda, cfg.max_pix_mcr);
fprintf('Subsample: %d pixels x %d wavelengths\n', size(D, 1), size(D, 2));
Dsimp = D(randperm(size(D, 1), min(cfg.n_simplisma, size(D, 1))), :);

sv = svd(D, 'econ');
var_sv = 100 * sv.^2 / sum(sv.^2);
figure('Color', 'w');
subplot(1, 2, 1); semilogy(1:15, sv(1:15), 'ko-', 'MarkerFaceColor', 'k');
xlabel('Component'); ylabel('Singular value (log scale)'); grid on;
subplot(1, 2, 2); plot(1:15, cumsum(var_sv(1:15)), 'ko-', 'MarkerFaceColor', 'k');
xlabel('Number of components'); ylabel('Cumulative explained variance (%)'); grid on;
exportgraphics(gcf, fullfile(cfg.output_folder, 'FigS1_singular_values.png'), 'Resolution', 300);
close(gcf);

ncs = 3:5;
res_nc = table('Size', [numel(ncs) 6], ...
    'VariableTypes', {'double','double','double','double','logical','double'}, ...
    'VariableNames', {'nc','lack_of_fit_pct','r2_pct','n_iter','converged','r2_svd_pct'});
figure('Color', 'w', 'Position', [100 100 1200 350]);
for i = 1:numel(ncs)
    nc = ncs(i);
    % SIMPLISMA + 5 random starts (and the full-model profiles for nc = 4);
    % the start with the lowest lack of fit is kept
    [~, S, info] = mcr_als_nn(D, simplisma_spectra(Dsimp, nc, cfg.simplisma_offset), cfg.als);
    for rep = 1:5
        [~, Sr, infor] = mcr_als_nn(D, Dsimp(randperm(size(Dsimp, 1), nc), :), cfg.als);
        if infor.lof < info.lof, S = Sr; info = infor; end
    end
    if nc == cfg.nc
        [~, Sp, infop] = mcr_als_nn(D, soptTots, cfg.als);
        if infop.lof <= info.lof + 0.01, S = Sp; info = infop; end
        S = match_components(S, soptTots);
    end
    res_nc(i, :) = {nc, info.lof, info.r2, info.n_iter, info.converged, sum(var_sv(1:nc))};
    subplot(1, numel(ncs), i);
    plot(lambda_nm, S', 'LineWidth', 1.3);
    title(sprintf('%d components (LOF %.2f %%)', nc, info.lof));
    xlabel('Wavelength (nm)'); ylabel('Normalised intensity'); grid on;
    xlim([min(lambda_nm) max(lambda_nm)]);
    legend(compose('C%d', 1:nc), 'Location', 'best');
end
exportgraphics(gcf, fullfile(cfg.output_folder, 'FigS2_profiles_3_4_5_components.png'), 'Resolution', 300);
close(gcf);
disp(res_nc);

[~, S_chk, info_chk] = mcr_als_nn(D, soptTots, cfg.als);
r_chk = diag(corr(S_chk', soptTots'))';
fprintf('Started from the full-model profiles: r2 = %.2f %%, correlation = %s\n', info_chk.r2, mat2str(r_chk, 3));

r_ini = nan(n_rep, cfg.nc); r_sub = nan(n_rep, cfg.nc);
S_ini = nan(cfg.nc, numel(cfg.idx_lambda), n_rep);
for k = 1:n_rep
    [~, S] = mcr_als_nn(D, Dsimp(randperm(size(Dsimp, 1), cfg.nc), :), cfg.als);   % random initialisation
    [S, ~, r_ini(k, :)] = match_components(S, soptTots);
    S_ini(:, :, k) = S;
    [Dk, ~] = extract_pixels(result, T, 1:height(T), cfg.idx_lambda, cfg.max_pix_mcr);   % other subsample
    S0 = simplisma_spectra(Dk(randperm(size(Dk, 1), min(cfg.n_simplisma, size(Dk, 1))), :), ...
        cfg.nc, cfg.simplisma_offset);
    [~, S] = mcr_als_nn(Dk, S0, cfg.als);
    [~, ~, r_sub(k, :)] = match_components(S, soptTots);
end
figure('Color', 'w', 'Position', [100 100 1000 600]);
for c = 1:cfg.nc
    subplot(2, 2, c); hold on;
    Sc = squeeze(S_ini(c, :, :));
    fill([lambda_nm fliplr(lambda_nm)], [min(Sc, [], 2)' fliplr(max(Sc, [], 2)')], ...
        [0.8 0.8 0.8], 'EdgeColor', 'none');
    plot(lambda_nm, soptTots(c, :) / norm(soptTots(c, :)), 'k', 'LineWidth', 1.5);
    title(sprintf('Component %d (median r = %.3f)', c, median(r_ini(:, c))));
    xlabel('Wavelength (nm)'); ylabel('Normalised intensity'); grid on;
    xlim([min(lambda_nm) max(lambda_nm)]);
end
legend({sprintf('Range over %d random initialisations', n_rep), 'Full-model profile'}, 'Location', 'best');
exportgraphics(gcf, fullfile(cfg.output_folder, 'FigS3_initialisation_stability.png'), 'Resolution', 300);
close(gcf);

stab = table((1:cfg.nc)', median(r_ini)', min(r_ini)', median(r_sub)', min(r_sub)', ...
    'VariableNames', {'component','r_median_random_init','r_min_random_init', ...
                      'r_median_subsamples','r_min_subsamples'});
disp(stab);
fx = fullfile(cfg.output_folder, 'TableS_mcr.xlsx');
writetable(res_nc, fx, 'Sheet', 'n_components');
writetable(stab,   fx, 'Sheet', 'stability');
writetable(table(info_chk.r2, r_chk, 'VariableNames', {'r2_pct','corr_with_soptTots'}), fx, 'Sheet', 'full_model_check');
end

function step_leave_one_pot_out(cfg)
% Leave-one-pot-out validation. In each fold the left-out pot is excluded
% from all model-building steps: MCR-ALS refit (initial estimates
% cfg.mcr_init), PCA reference model on the mock-inoculated calibration
% leaves (all dates), control limit, and pixel-level model. Each leaf x date
% observation receives exactly one out-of-sample prediction.
load(cfg.data_file, 'result', 'soptTots');
load(fullfile(cfg.output_folder, 'leaf_table.mat'), 'T');
rng(cfg.seed);
run_mode = mode_name(cfg);
fprintf('\n===== Leave-one-pot-out validation (%s) =====\n', run_mode);

% Mean spectrum of each leaf x date (Eq. 6 is linear: mean(D)*pinv(S) = mean(D*pinv(S)))
nObs = height(T);
Dmean = zeros(nObs, numel(cfg.idx_lambda));
for r = 1:nObs
    Dmean(r, :) = mean(result{T.i_res(r)}.data(T.row_start(r):T.row_end(r), cfg.idx_lambda), 1);
end

pots = unique(T.pot);
nF = numel(pots);
Cm_oos = nan(nObs, cfg.nc);
T2 = nan(nObs, 1); Q = nan(nObs, 1); T2lim = nan(nObs, 1); Qlim = nan(nObs, 1); fold = nan(nObs, 1);
models = struct('pot_id', {}, 'S', {}, 'M', {}, 'Mpix', {}, 'info', {}, 'r_ref', {});
for f = 1:nF
    is_val = T.pot == pots(f);
    is_cal = ~is_val;
    % MCR-ALS on the calibration pots
    if cfg.mcr_mock_only, rows_mcr = find(is_cal & ~T.inoculated); else, rows_mcr = find(is_cal); end
    [Dcal, id_cal] = extract_pixels(result, T, rows_mcr, cfg.idx_lambda, cfg.max_pix_mcr);
    if strcmp(cfg.mcr_init, 'full_model')
        S0 = soptTots;
    else
        Dsimp = Dcal(randperm(size(Dcal, 1), min(cfg.n_simplisma, size(Dcal, 1))), :);
        S0 = simplisma_spectra(Dsimp, cfg.nc, cfg.simplisma_offset);
    end
    [~, S, info] = mcr_als_nn(Dcal, S0, cfg.als);
    [S, ~, r_ref] = match_components(S, soptTots);     % same component order as the full model
    % Mean contributions of all leaves with the profiles of the fold
    if strcmp(cfg.projection, 'pinv')
        Cm = Dmean * pinv(S);
    else
        Cm = zeros(nObs, cfg.nc);                       % NNLS is not linear: per pixel, then mean
        for r = 1:nObs
            Dp = result{T.i_res(r)}.data(T.row_start(r):T.row_end(r), cfg.idx_lambda);
            Cm(r, :) = mean(project_spectra(Dp, S, 'nnls'), 1);
        end
    end
    % PCA reference model: mock-inoculated calibration leaves
    ref = is_cal & ~T.inoculated;
    if cfg.mspc.ref_by_genotype
        ref = ref & strcmp(T.genotype, T.genotype{find(is_val, 1)});
    end
    M = mspc_fit(Cm(ref, :), cfg.mspc);
    % Left-out pot
    [t2v, qv] = mspc_apply(M, Cm(is_val, :));
    Cm_oos(is_val, :) = Cm(is_val, :);
    T2(is_val) = t2v; Q(is_val) = qv;
    switch cfg.mspc.limit
        case 'F',    T2lim(is_val) = M.T2lim_F;
        case 'chi2', T2lim(is_val) = M.T2lim_chi2;
        case 'emp',  T2lim(is_val) = M.T2lim_emp;
    end
    Qlim(is_val) = M.Qlim;
    fold(is_val) = f;
    % Pixel-level model (maps): pixels of mock-inoculated calibration leaves,
    % empirical limit (percentile cfg.pix.percentile)
    pix_ref = ref(id_cal);
    Cpix = project_spectra(Dcal(pix_ref, :), S, cfg.projection);
    Mpix = mspc_fit(Cpix, struct('k', M.k, 'alpha', cfg.mspc.alpha));
    Mpix.T2lim_pix = prctile(mspc_apply(Mpix, Cpix), cfg.pix.percentile);
    models(f) = struct('pot_id', pots(f), 'S', S, 'M', M, 'Mpix', Mpix, 'info', info, 'r_ref', r_ref);
    fprintf('Fold %2d/%d (pot %d): r2 = %.2f %%, k = %d, corr. with full model = %s\n', ...
        f, nF, pots(f), info.r2, M.k, mat2str(r_ref, 3));
end

R = T(:, {'pot','date','dpi','leaf','leaf_id','image_id','genotype','inoculated'});
for c = 1:cfg.nc, R.(sprintf('C%d', c)) = Cm_oos(:, c); end
R.T2 = T2; R.T2lim = T2lim; R.Q = Q; R.Qlim = Qlim; R.fold = fold;
R.alarm_T2 = R.T2 > R.T2lim;
if cfg.mspc.use_Q && all(~isnan(R.Qlim))
    R.alarm = R.alarm_T2 | R.Q > R.Qlim;
else
    R.alarm = R.alarm_T2;                               % with k = nc PCs, SPE = 0
end

k_f = arrayfun(@(m) m.M.k, models); N_f = arrayfun(@(m) m.M.N, models);
lim_f = arrayfun(@(m) m.M.T2lim_F, models); r2_f = arrayfun(@(m) m.info.r2, models);
expl_f = cell2mat(arrayfun(@(m) m.M.expl(1:max(k_f))', models(:), 'UniformOutput', false));
fprintf('MCR-ALS r2 %.2f-%.2f %%; PCs %s; reference observations %d-%d; T2 limit %.2f-%.2f\n', ...
    min(r2_f), max(r2_f), mat2str(unique(k_f)), min(N_f), max(N_f), min(lim_f), max(lim_f));
fprintf('Variance explained per PC: first fold %s %%; range over folds %s-%s %%\n', ...
    mat2str(models(1).M.expl(1:max(k_f))', 3), mat2str(min(expl_f, [], 1), 3), mat2str(max(expl_f, [], 1), 3));
fprintf('Minimum correlation of the fold profiles with the full model: %s\n', ...
    mat2str(min(vertcat(models.r_ref), [], 1), 3));

[Met, boot_dates] = detection_metrics(R, cfg.n_boot, cfg.seed, cfg.mspc.alpha);
[Det, Tests] = detection_table(R, boot_dates, cfg.mspc.alpha);
disp(Det);
mock = ~R.inoculated;
fprintf('False alarms in mock-inoculated leaves: %d of %d (%.1f %%); 10-14 dpi %d of %d; 16-17 dpi %d of %d\n', ...
    sum(R.alarm(mock)), sum(mock), 100 * mean(R.alarm(mock)), ...
    sum(R.alarm(mock & R.dpi >= 10 & R.dpi <= 14)), sum(mock & R.dpi >= 10 & R.dpi <= 14), ...
    sum(R.alarm(mock & R.dpi >= 16)), sum(mock & R.dpi >= 16));
[gm, pot_k, gen] = findgroups(R.pot(mock), R.genotype(mock));
FA_pot = table(pot_k, gen, splitapply(@numel, R.alarm(mock), gm), ...
    round(100 * splitapply(@mean, R.alarm(mock), gm), 1), ...
    round(100 * splitapply(@mean, R.alarm(mock & R.dpi <= 0), findgroups(R.pot(mock & R.dpi <= 0))), 1), ...
    'VariableNames', {'pot_id','genotype','n_obs','pct_false_alarms','pct_before_inoculation'});
FA_pot = sortrows(FA_pot, 'pct_false_alarms', 'descend');
disp(FA_pot);
Stab = array2table([[models.pot_id]' vertcat(models.r_ref)], ...
    'VariableNames', [{'validation_pot'}, compose('r_C%d', 1:cfg.nc)]);

fx = fullfile(cfg.output_folder, sprintf('Table_metrics_%s.xlsx', run_mode));
writetable(Met,    fx, 'Sheet', 'metrics_genotype_dpi');
writetable(Det,    fx, 'Sheet', 'detection_date');
writetable(Tests,  fx, 'Sheet', 'fisher_tests_per_dpi');
writetable(R,      fx, 'Sheet', 'observations');
writetable(FA_pot, fx, 'Sheet', 'false_alarms_per_pot');
writetable(Stab,   fx, 'Sheet', 'profile_stability');
save(fullfile(cfg.output_folder, sprintf('validation_%s.mat', run_mode)), ...
    'R', 'models', 'Met', 'Det', 'Tests', 'Dmean', 'cfg');
plot_metrics(R, Met, cfg, fullfile(cfg.output_folder, sprintf('FigS_metrics_%s.png', run_mode)), '');
end

function [Det, Tests] = detection_table(R, boot_dates, alpha)
genos = unique(R.genotype);
Det = table(); Tests = table();
for g = 1:numel(genos)
    [d0, tp] = detection_date(R, genos{g}, alpha);
    ci = prctile(boot_dates(:, g), [2.5 50 97.5]);
    Det = [Det; table(genos(g), d0, ci(2), ci(1), ci(3), mean(isnan(boot_dates(:, g))), ...
        'VariableNames', {'genotype','detection_dpi','boot_median','CI_low','CI_high','frac_boot_no_detection'})]; %#ok<AGROW>
    tp.genotype = repmat(genos(g), height(tp), 1);
    Tests = [Tests; tp]; %#ok<AGROW>
end
end

function plot_metrics(R, Met, cfg, file, ttl)
% Sensitivity and false-positive rate per dpi (check figure, not in the article)
col_inoc = [0.80 0.15 0.15]; col_mock = [0.25 0.25 0.25];
genos = unique(R.genotype);
figure('Color', 'w', 'Position', [100 100 1100 420]);
for g = 1:numel(genos)
    subplot(1, numel(genos), g); hold on;
    m = Met(strcmp(Met.genotype, genos{g}), :);
    errorbar(m.dpi, m.sensitivity, m.sensitivity - m.sens_CI_low, m.sens_CI_high - m.sensitivity, ...
        'o-', 'Color', col_inoc, 'MarkerFaceColor', col_inoc);
    errorbar(m.dpi, m.FPR, m.FPR - m.FPR_CI_low, m.FPR_CI_high - m.FPR, ...
        's-', 'Color', col_mock, 'MarkerFaceColor', col_mock);
    ylim([0 1]); grid on; box on;
    xlabel('Days post-inoculation'); ylabel('Proportion of leaves');
    xlim([min(cfg.dpi) - 0.6, max(cfg.dpi) + 0.6]); xticks([-1 0 4 6 8 10 12 14 16]);
    title(strtrim([genos{g} ' ' ttl]));
end
legend({'Sensitivity', 'False-positive rate'}, 'Location', 'northwest');
exportgraphics(gcf, file, 'Resolution', 300);
close(gcf);
end

function step_mixed_models(cfg)
% Linear mixed models on the out-of-sample contributions and log T2:
% y ~ genotype * treatment * dpi + (1|pot) + (1|leaf), and, at each dpi and
% genotype, y ~ treatment + (1|pot) with Holm correction (0-17 dpi);
% sustained divergence = first dpi from which the difference is significant
% with the same sign at all later dpi.
run_mode = mode_name(cfg);
fprintf('\n===== Mixed models (%s) =====\n', run_mode);
load(fullfile(cfg.output_folder, sprintf('validation_%s.mat', run_mode)), 'R');
R.treat  = categorical(double(R.inoculated), [0 1], {'mock-inoculated','inoculated'});
R.geno   = categorical(R.genotype);
R.dpi_c  = categorical(R.dpi);
R.pot_id = categorical(R.pot);
R.leaf_c = categorical(R.leaf_id);
R.logT2  = log(R.T2);
genos = categories(R.geno);
dpis  = unique(R.dpi);
vars  = [compose('C%d', 1:cfg.nc), {'logT2'}];

Anova_all = table();
for v = 1:numel(vars)
    R.y = R.(vars{v});
    lme = fitlme(R, 'y ~ geno*treat*dpi_c + (1|pot_id) + (1|leaf_c)');
    a = anova(lme);
    if isa(a, 'dataset'), A = dataset2table(a); else, A = a; end
    A.variable = repmat(vars(v), height(A), 1);
    Anova_all = [Anova_all; A]; %#ok<AGROW>
end
disp(Anova_all(contains(Anova_all.Term, 'treat:geno'), :));

Effects = table(); Diverg = table();
for v = 1:numel(vars)
    for g = 1:numel(genos)
        est = nan(numel(dpis), 1); lo = est; hi = est; p = est;
        for d = 1:numel(dpis)
            sub = R(R.geno == genos{g} & R.dpi == dpis(d), :);
            if numel(unique(sub.treat)) < 2, continue; end
            sub.y = sub.(vars{v});
            try
                lme = fitlme(sub, 'y ~ treat + (1|pot_id)');
            catch ME
                warning('%s, %s, %d dpi: %s', vars{v}, genos{g}, dpis(d), ME.message);
                continue;
            end
            row = strcmp(lme.Coefficients.Name, 'treat_inoculated');
            ci = coefCI(lme);
            est(d) = lme.Coefficients.Estimate(row);
            lo(d) = ci(row, 1); hi(d) = ci(row, 2);
            p(d) = lme.Coefficients.pValue(row);
        end
        padj = nan(size(p));
        post = dpis >= 0;
        padj(post) = holm(p(post));
        Effects = [Effects; table(repmat(vars(v), numel(dpis), 1), repmat(genos(g), numel(dpis), 1), ...
            dpis, est, lo, hi, p, padj, 'VariableNames', ...
            {'variable','genotype','dpi','diff_inoc_minus_mock','CI_low','CI_high','p','p_holm'})]; %#ok<AGROW>
        k = find(post & padj < cfg.mspc.alpha, 1);
        if isempty(k), dd = NaN; else, dd = dpis(k); end
        ipost = find(post);
        sig = padj(ipost) < cfg.mspc.alpha;
        sg  = sign(est(ipost));
        ds = NaN; sgn_val = NaN;
        for j = 1:numel(ipost)
            if all(sig(j:end)) && all(sg(j:end) == sg(j))
                ds = dpis(ipost(j)); sgn_val = sg(j);
                break;
            end
        end
        Diverg = [Diverg; table(vars(v), genos(g), dd, ds, sgn_val, 'VariableNames', ...
            {'variable','genotype','first_signif_dpi','sustained_divergence_dpi','sgn_val'})]; %#ok<AGROW>
    end
end
disp(Diverg);
fx = fullfile(cfg.output_folder, sprintf('Table_mixed_models_%s.xlsx', run_mode));
writetable(Anova_all, fx, 'Sheet', 'anova');
writetable(Effects,   fx, 'Sheet', 'treatment_effect_per_dpi');
writetable(Diverg,    fx, 'Sheet', 'divergence_dpi');
end

function step_pixels(cfg)
% Percentage of pixels per leaf assigned to the infection-associated class:
% pixels projected onto the profiles of the fold in which their pot was left
% out, T2 in the pixel-level model of that fold, assigned when above the
% empirical limit (99th percentile of the reference pixels).
run_mode = mode_name(cfg);
fprintf('\n===== Pixel assignment (%s) =====\n', run_mode);
load(cfg.data_file, 'result');
load(fullfile(cfg.output_folder, 'leaf_table.mat'), 'T');
load(fullfile(cfg.output_folder, sprintf('validation_%s.mat', run_mode)), 'models');
nObs = height(T);
pct = nan(nObs, 1);
T2esc = cell(nObs, 1);
pot_mod = [models.pot_id];
for r = 1:nObs
    Dp = result{T.i_res(r)}.data(T.row_start(r):T.row_end(r), cfg.idx_lambda);
    md = models(pot_mod == T.pot(r));
    t2b = mspc_apply(md.Mpix, project_spectra(Dp, md.S, cfg.projection));
    pct(r) = 100 * mean(t2b > md.Mpix.T2lim_pix);
    T2esc{r} = t2b / md.Mpix.T2lim_pix;
end
P = T(:, {'pot','date','dpi','leaf','leaf_id','genotype','inoculated','npix'});
P.pct_assigned = pct;
[g, geno, inoc, dpi] = findgroups(P.genotype, P.inoculated, P.dpi);
summ = table(geno, inoc, dpi, splitapply(@numel, pct, g), splitapply(@median, pct, g), ...
    splitapply(@(x) prctile(x, 75), pct, g), splitapply(@max, pct, g), splitapply(@(x) mean(x > 20), pct, g), ...
    'VariableNames', {'genotype','inoculated','dpi','n_leaves','median_pct','p75_pct','max_pct','frac_leaves_over_20pct'});
summ = sortrows(summ, {'genotype','inoculated','dpi'});
mock = ~P.inoculated;
fprintf('Mock-inoculated leaves: median %.1f %%, maximum %.1f %%; > 20 %% in %d of %d observations (%d of %d at 16-17 dpi)\n', ...
    median(pct(mock)), max(pct(mock)), sum(pct(mock) > 20), sum(mock), ...
    sum(pct(mock & P.dpi >= 16) > 20), sum(mock & P.dpi >= 16));
fprintf('Inoculated leaves before inoculation: maximum %.1f %%\n', max(pct(~mock & P.dpi <= 0)));
fx = fullfile(cfg.output_folder, sprintf('Table_pixels_%s.xlsx', run_mode));
writetable(P,    fx, 'Sheet', 'per_leaf_date');
writetable(summ, fx, 'Sheet', 'summary_group_dpi');
info_table_pix = table(T.pot, T.date, T.leaf, T2esc, cellfun(@mean, T2esc), pct, ...
    'VariableNames', {'pot','date','leaf','T2_all','T2_mean','pct_bad_pixels'});
save(fullfile(cfg.output_folder, sprintf('info_table_pix_%s.mat', run_mode)), 'info_table_pix');
end

function step_time_updated_reference(cfg)
% Reference built from mock-inoculated calibration leaves of the same date
% ('concurrent') or of the two preceding dates ('preceding'); everything
% else as in the main validation (fold profiles, 4 PCs, F limit).
run_mode = mode_name(cfg);
assert(strcmp(cfg.projection, 'pinv'), 'The time-updated references use the linear projection.');
load(fullfile(cfg.output_folder, sprintf('validation_%s.mat', run_mode)), 'models', 'Dmean');
load(fullfile(cfg.output_folder, 'leaf_table.mat'), 'T');
alpha = cfg.mspc.alpha;
dpis = cfg.dpi;
for variant = {'preceding', 'concurrent'}
    fprintf('\n===== Time-updated reference: %s =====\n', variant{1});
    nObs = height(T);
    T2 = nan(nObs, 1); T2lim = nan(nObs, 1); Nref = nan(nObs, 1); Cm_oos = nan(nObs, cfg.nc);
    pot_mod = [models.pot_id];
    for f = 1:numel(models)
        is_val = T.pot == pot_mod(f);
        is_cal = ~is_val;
        Cm = Dmean * pinv(models(f).S);
        Cm_oos(is_val, :) = Cm(is_val, :);
        for j = 1:numel(dpis)
            if strcmp(variant{1}, 'preceding')
                if j == 1, ref_dates = dpis(1); else, ref_dates = dpis(max(1, j - 2):j - 1); end
            else
                ref_dates = dpis(j);
            end
            ref = is_cal & ~T.inoculated & ismember(T.dpi, ref_dates);
            val = is_val & T.dpi == dpis(j);
            if ~any(val), continue; end
            M = mspc_fit(Cm(ref, :), struct('k', 4, 'alpha', alpha));
            T2(val) = mspc_apply(M, Cm(val, :));
            T2lim(val) = M.T2lim_F;
            Nref(val) = M.N;
        end
    end
    R = T(:, {'pot','date','dpi','leaf','leaf_id','image_id','genotype','inoculated'});
    for c = 1:cfg.nc, R.(sprintf('C%d', c)) = Cm_oos(:, c); end
    R.T2 = T2; R.T2lim = T2lim; R.Nref = Nref;
    R.alarm_T2 = R.T2 > R.T2lim;
    R.alarm = R.alarm_T2;
    fprintf('Reference observations per model: %d-%d\n', min(Nref), max(Nref));
    [Met, boot_dates] = detection_metrics(R, cfg.n_boot, cfg.seed, alpha);
    [Det, Tests] = detection_table(R, boot_dates, alpha);
    disp(Det);
    mock = ~R.inoculated;
    fprintf('False alarms: %.1f %% overall, %.1f %% at 10-14 dpi, %.1f %% at 16-17 dpi\n', ...
        100 * mean(R.alarm(mock)), 100 * mean(R.alarm(mock & R.dpi >= 10 & R.dpi <= 14)), ...
        100 * mean(R.alarm(mock & R.dpi >= 16)));
    suffix = sprintf('%s_mw_%s', run_mode, variant{1});
    fx = fullfile(cfg.output_folder, sprintf('Table_metrics_%s.xlsx', suffix));
    writetable(Met,   fx, 'Sheet', 'metrics_genotype_dpi');
    writetable(Det,   fx, 'Sheet', 'detection_date');
    writetable(Tests, fx, 'Sheet', 'fisher_tests_per_dpi');
    writetable(R,     fx, 'Sheet', 'observations');
    save(fullfile(cfg.output_folder, sprintf('validation_%s.mat', suffix)), 'R', 'Met', 'Det');
    plot_metrics(R, Met, cfg, fullfile(cfg.output_folder, sprintf('FigS_metrics_%s.png', suffix)), ...
        sprintf('(%s reference)', variant{1}));
end
end

function step_export_metadata(cfg, excluded)
% metadata/leaves_by_date.csv, folds.csv and excluded_pots.csv
run_mode = mode_name(cfg);
load(fullfile(cfg.output_folder, 'leaf_table.mat'), 'T');
load(fullfile(cfg.output_folder, sprintf('validation_%s.mat', run_mode)), 'R');
dir_meta = fullfile(cfg.output_folder, 'metadata');
if ~exist(dir_meta, 'dir'), mkdir(dir_meta); end
treatment = repmat({'mock'}, height(T), 1);
treatment(T.inoculated) = {'inoculated'};
L = table(T.pot, T.side, T.leaf, T.leaf_id, T.date, T.dpi, T.genotype, treatment, T.npix, T.image_id, R.fold, ...
    'VariableNames', {'pot','side','leaf','leaf_id','date','dpi','genotype','treatment','n_pixels','image_id','fold'});
writetable(L, fullfile(dir_meta, 'leaves_by_date.csv'));
[pots, ia] = unique(L.pot);
writetable(table(pots, L.genotype(ia), L.treatment(ia), L.fold(ia), ...
    'VariableNames', {'pot','genotype','treatment','fold'}), fullfile(dir_meta, 'folds.csv'));
reason = arrayfun(@(e) sprintf(['Leaves outside the healthy reference before inoculation ' ...
    '(one-sided binomial test at -1 and 0 dpi; iteration %d of the check)'], e.iteration), excluded(:), 'UniformOutput', false);
writetable(table([excluded.pot]', {excluded.genotype}', repmat({'mock'}, numel(excluded), 1), reason, ...
    'VariableNames', {'pot','genotype','treatment','reason'}), fullfile(dir_meta, 'excluded_pots.csv'));
copyfile(fullfile(cfg.output_folder, 'excluded_leaf_observations.csv'), fullfile(dir_meta, 'excluded_leaf_observations.csv'));
copyfile(fullfile(cfg.output_folder, 'excluded_leaves_low_confidence.csv'), fullfile(dir_meta, 'excluded_leaves_low_confidence.csv'));
fprintf('\nMetadata: %d observations, %d leaves, %d pots (%s)\n', height(L), numel(unique(L.leaf_id)), numel(pots), dir_meta);
end


%% ========================================================================
%  STEP 'sensitivity': summary (Table S3)
%  ========================================================================
function step_sensitivity_summary(cfg)
% Compares each sensitivity analysis with the main analysis: detection
% dates, false alarms, agreement of the alarm decisions and of log T2 (same
% observations), and correlation of the fold profiles.
main_mode = mode_name(cfg);
main = load(fullfile(cfg.output_folder, sprintf('validation_%s.mat', main_mode)), 'R', 'Det', 'models');
c3 = cfg; c3.mspc.k = [];
cs = cfg; cs.mcr_init = 'simplisma';
cn = cfg; cn.projection = 'nnls';
runs = {'Main analysis (four PCs)',       fullfile(cfg.output_folder, sprintf('validation_%s.mat', main_mode));
        'Three PCs, T2 and SPE',          fullfile(cfg.output_folder, sprintf('validation_%s.mat', mode_name(c3)));
        'SIMPLISMA initial estimates in each fold', fullfile(cfg.output_folder, sprintf('validation_%s.mat', mode_name(cs)));
        'Non-negative projection',        fullfile(cfg.output_folder, sprintf('validation_%s.mat', mode_name(cn)));
        'All 20 pots',                    fullfile(cfg.output_folder, 'all_20_pots', sprintf('validation_%s.mat', main_mode));
        'Reference of the same date',     fullfile(cfg.output_folder, sprintf('validation_%s_mw_concurrent.mat', main_mode));
        'Reference of the two preceding dates', fullfile(cfg.output_folder, sprintf('validation_%s_mw_preceding.mat', main_mode))};
S = table();
for i = 1:size(runs, 1)
    if exist(runs{i, 2}, 'file') ~= 2
        fprintf('Not found (skipped): %s\n', runs{i, 2});
        continue
    end
    v = load(runs{i, 2});
    R = v.R;
    mock = ~R.inoculated;
    d67 = v.Det.detection_dpi(strcmp(v.Det.genotype, 'EPO67'));
    d68 = v.Det.detection_dpi(strcmp(v.Det.genotype, 'EPO68'));
    % agreement with the main analysis on the common observations
    [~, ia, ib] = intersect(strcat(R.leaf_id, '_', R.date), strcat(main.R.leaf_id, '_', main.R.date), 'stable');
    same_alarm = 100 * mean(R.alarm(ia) == main.R.alarm(ib));
    r_logT2 = corr(log(R.T2(ia)), log(main.R.T2(ib)));
    if isfield(v, 'models') && numel(v.models) == numel(main.models) && isequal([v.models.pot_id], [main.models.pot_id])
        rS = nan(numel(v.models), cfg.nc);
        for f = 1:numel(v.models)
            rS(f, :) = diag(corr(v.models(f).S', main.models(f).S'))';
        end
        min_r_profiles = mat2str(min(rS, [], 1), 2);
    else
        min_r_profiles = '';
    end
    S = [S; table(runs(i, 1), d67, d68, 100 * mean(R.alarm(mock)), ...
        100 * mean(R.alarm(mock & R.dpi >= 10 & R.dpi <= 14)), 100 * mean(R.alarm(mock & R.dpi >= 16)), ...
        same_alarm, r_logT2, {min_r_profiles}, 'VariableNames', ...
        {'analysis','detection_dpi_EPO67','detection_dpi_EPO68','pct_false_alarms', ...
         'pct_false_alarms_10_14dpi','pct_false_alarms_16_17dpi','pct_same_alarm_as_main', ...
         'r_logT2_with_main','min_r_fold_profiles_with_main'})]; %#ok<AGROW>
end
disp(S);
writetable(S, fullfile(cfg.output_folder, 'TableS3_S4_sensitivity_summary.xlsx'));
end


%% ========================================================================
%  STEP 'figures'
%  ========================================================================
function fig2_trajectories(cfg)
% Figure 2: profile of each component (mean and range over folds) and mean
% contribution per group and dpi, with 95 % CI on pot means. Full page.
run_mode = mode_name(cfg);
load(fullfile(cfg.output_folder, sprintf('validation_%s.mat', run_mode)), 'R', 'models');
R.treat = categorical(double(R.inoculated), [0 1], {'mock-inoculated','inoculated'});
R.geno  = categorical(R.genotype);
dpis = unique(R.dpi);
lambda_nm = cfg.lambda_nm;
labels = {'Component 1 (healthy-leaf-like profile)', 'Component 2 (senescence-associated profile)', ...
          'Component 3 (visible-absorption profile)', 'Component 4 (infection-associated profile)'};
groups = {'EPO67','mock-inoculated';'EPO67','inoculated';'EPO68','mock-inoculated';'EPO68','inoculated'};
cols = [0.55 0.55 0.55; 0.85 0.20 0.20; 0.15 0.15 0.15; 0.95 0.60 0.60];
markers = {'s','o','s','o'};
S_all = cat(3, models.S);
figure('Color', 'w', 'Position', [50 50 900 1180]);
for c = 1:cfg.nc
    subplot(cfg.nc, 2, 2*c - 1); hold on;
    Sc = squeeze(S_all(c, :, :));
    fill([lambda_nm fliplr(lambda_nm)], [min(Sc, [], 2)' fliplr(max(Sc, [], 2)')], [0.85 0.85 0.85], 'EdgeColor', 'none');
    plot(lambda_nm, mean(Sc, 2), 'k', 'LineWidth', 1.4);
    title(labels{c}, 'FontSize', 10);
    if c == cfg.nc, xlabel('Wavelength (nm)'); end
    ylabel('Normalised intensity'); grid on; box on;
    xlim([min(lambda_nm) max(lambda_nm)]);

    subplot(cfg.nc, 2, 2*c); hold on;
    y = R.(sprintf('C%d', c));
    legend_txt = {};
    for k = 1:size(groups, 1)
        s = R.geno == groups{k, 1} & R.treat == groups{k, 2};
        m = nan(numel(dpis), 1); h = m; npot = 0;
        for d = 1:numel(dpis)
            sd = s & R.dpi == dpis(d);
            pot_k = unique(R.pot(sd));
            if isempty(pot_k), continue; end
            means = arrayfun(@(q) mean(y(sd & R.pot == q)), pot_k);
            m(d) = mean(means);
            if numel(means) > 1, h(d) = tinv(0.975, numel(means) - 1) * std(means) / sqrt(numel(means)); end
            npot = max(npot, numel(pot_k));
        end
        errorbar(dpis + (k - 2.5) * 0.12, m, h, ['-' markers{k}], 'Color', cols(k, :), ...
            'MarkerFaceColor', cols(k, :), 'MarkerSize', 4, 'CapSize', 3);
        legend_txt{end+1} = sprintf('%s %s (%d pots)', groups{k, 1}, groups{k, 2}, npot); %#ok<AGROW>
    end
    if c == cfg.nc, xlabel('Days post-inoculation'); end
    ylabel('Mean contribution'); grid on; box on;
    xlim([min(cfg.dpi) - 0.6, max(cfg.dpi) + 0.6]); xticks([0 4 8 12 16]);
    if c == 2, legend(legend_txt, 'Location', 'northwest', 'FontSize', 8); end   % empty area of this panel
end
exportgraphics(gcf, fullfile(cfg.output_folder, sprintf('Fig2_trajectories_%s.png', run_mode)), 'Resolution', 300);
close(gcf);
end

function fig3_loadings(cfg)
% Figure 3: loadings of the PCA reference model. Bars: first fold; error
% bars: range over all folds (sign of each PC aligned with the first fold).
run_mode = mode_name(cfg);
load(fullfile(cfg.output_folder, sprintf('validation_%s.mat', run_mode)), 'models');
nf = numel(models);
P1 = models(1).M.P;
[nc, k] = size(P1);
% sign convention (the sign of a PC is arbitrary): largest loading positive
for j = 1:k
    [~, im] = max(abs(P1(:, j)));
    if P1(im, j) < 0, P1(:, j) = -P1(:, j); end
end
% the sign of a PC is arbitrary: each PC is oriented so that its largest
% loading (in absolute value) is positive
[~, im] = max(abs(P1), [], 1);
P1 = P1 .* sign(P1(sub2ind(size(P1), im, 1:k)));
Pall = nan(nc, k, nf); expl = nan(k, nf); r = nan(k, nf);
for f = 1:nf
    Pf = models(f).M.P;
    for j = 1:k
        s = sign(Pf(:, j)' * P1(:, j)); if s == 0, s = 1; end
        Pall(:, j, f) = s * Pf(:, j);
        r(j, f) = abs(corr(Pf(:, j), P1(:, j)));
    end
    expl(:, f) = models(f).M.expl(1:k);
end
Pmin = min(Pall, [], 3); Pmax = max(Pall, [], 3);
figure('Color', 'w', 'Position', [100 100 800 380]);
b = bar(P1, 'grouped', 'EdgeColor', 'none');
greys = linspace(0.15, 0.80, k);
hold on;
for j = 1:k
    b(j).FaceColor = greys(j) * [1 1 1];
    errorbar(b(j).XEndPoints, P1(:, j), P1(:, j) - Pmin(:, j), Pmax(:, j) - P1(:, j), ...
        'k', 'LineStyle', 'none', 'CapSize', 3, 'HandleVisibility', 'off');
end
hold off;
set(gca, 'XTickLabel', compose('Component %d', 1:nc), 'FontSize', 11);
ylabel('Loading'); grid on; box on;
legend(b, arrayfun(@(j) sprintf('PC%d (%.1f%%)', j, models(1).M.expl(j)), 1:k, 'UniformOutput', false), ...
    'Location', 'eastoutside');
exportgraphics(gcf, fullfile(cfg.output_folder, sprintf('Fig3_loadings_%s.png', run_mode)), 'Resolution', 300);
close(gcf);
for j = 1:k
    fprintf('PC%d: variance %.1f-%.1f %% over %d folds; minimum |r| with the first fold %.3f\n', ...
        j, min(expl(j, :)), max(expl(j, :)), nf, min(r(j, :)));
end
end

function fig4_T2(cfg)
% Figure 4: out-of-sample T2 of each leaf (points) and median per group and
% dpi (lines), with the control limit.
run_mode = mode_name(cfg);
load(fullfile(cfg.output_folder, sprintf('validation_%s.mat', run_mode)), 'R');
rng(cfg.seed);
genos = unique(R.genotype);
col_inoc = [0.80 0.15 0.15]; col_mock = [0.25 0.25 0.25];
off = 0.18;
figure('Color', 'w', 'Position', [100 100 1100 460]);
tl = tiledlayout(1, numel(genos), 'TileSpacing', 'compact', 'Padding', 'compact');
hmed = gobjects(1, 2);
for g = 1:numel(genos)
    ax = nexttile(tl); hold on;
    for inoc = [false true]
        s = strcmp(R.genotype, genos{g}) & R.inoculated == inoc;
        if inoc, col = col_inoc; sgn = 1; else, col = col_mock; sgn = -1; end
        x = R.dpi(s) + sgn * off + (rand(sum(s), 1) - 0.5) * 0.15;
        scatter(x, R.T2(s), 10, col, 'filled', 'MarkerFaceAlpha', 0.30, 'HandleVisibility', 'off');
        dd = unique(R.dpi(s));
        med = arrayfun(@(d) median(R.T2(s & R.dpi == d)), dd);
        hmed(1 + inoc) = plot(dd + sgn * off, med, '-o', 'Color', col, 'MarkerFaceColor', col, ...
            'MarkerSize', 4, 'LineWidth', 1.8);
    end
    hlim = yline(mean(R.T2lim), '--k', 'LineWidth', 1);   % fold limits differ only in the second decimal
    set(ax, 'YScale', 'log', 'FontSize', 11);
    xlabel('Days post-inoculation'); ylabel('Hotelling''s T^2');
    xlim([min(cfg.dpi) - 0.6, max(cfg.dpi) + 0.6]); xticks([-1 0 4 6 8 10 12 14 16]);
    title(genos{g}); grid on; box on;
end
lg = legend([hmed hlim], {'Mock-inoculated (median)', 'Inoculated (median)', ...
    sprintf('Control limit (%.2f-%.2f)', min(R.T2lim), max(R.T2lim))}, 'Orientation', 'horizontal', 'FontSize', 11);
if isprop(lg, 'Layout')
    lg.Layout.Tile = 'south';            % below both panels (R2020b or later)
else
    lg.Location = 'southoutside';        % R2020a
end
exportgraphics(gcf, fullfile(cfg.output_folder, sprintf('Fig4_T2_validation_%s.png', run_mode)), 'Resolution', 300);
close(gcf);
end

function fig5_pixels(cfg)
% Figure 5: percentage of pixels per leaf assigned to the infection-associated class
run_mode = mode_name(cfg);
P = readtable(fullfile(cfg.output_folder, sprintf('Table_pixels_%s.xlsx', run_mode)), 'Sheet', 'per_leaf_date');
if ~islogical(P.inoculated), P.inoculated = logical(P.inoculated); end
genos = unique(P.genotype);
figure('Color', 'w', 'Position', [50 50 1300 380]);
for gg = 1:numel(genos)
    subplot(1, numel(genos), gg);
    s = strcmp(P.genotype, genos{gg});
    group = categorical(double(P.inoculated(s)), [0 1], {'Mock-inoculated','Inoculated'});
    boxchart(categorical(P.dpi(s)), P.pct_assigned(s), 'GroupByColor', group, 'MarkerStyle', '.', 'BoxWidth', 0.7);
    ylim([0 100]); grid on;
    xlabel('Days post-inoculation');
    ylabel({'Pixels assigned to the', 'infection-associated class (%)'});
    title(genos{gg});
    if gg == 1, legend('Location', 'northwest'); end
end
colororder([0.4 0.4 0.4; 0.85 0.2 0.2]);
exportgraphics(gcf, fullfile(cfg.output_folder, sprintf('Fig5_pixels_%s.png', run_mode)), 'Resolution', 300);
close(gcf);
end

function fig6_maps(cfg)
% Figure 6: (a) pixel maps (red = assigned to the infection-associated
% class, green = other leaf pixels); (b) colour composites of the same leaf
% pixels (bands near 640, 550 and 460 nm), cropped, padded and stretched
% exactly as the maps
run_mode = mode_name(cfg);
load(fullfile(cfg.output_folder, sprintf('info_table_pix_%s.mat', run_mode)), 'info_table_pix');
load(cfg.alldata_file, 'metadataG','metadataD','G_pixin1','G_pixin2','G_pixin3','G_pixout', ...
     'D_pixin1','D_pixin2','D_pixin3','D_pixout','G_Z1','G_Z2','G_Z3','D_Z1','D_Z2','D_Z3');
[~, wl] = full_model_profiles;
rgb_nm = [640 550 460];                  % red, green and blue channels
band_rgb = cell(1, 3);
for k = 1:3
    [~, ib] = sort(abs(wl - rgb_nm(k)));
    band_rgb{k} = sort(ib(1:11));        % mean of the 11 closest bands
end
rows = {32, 'G', sprintf('EPO67\ninoculated');
        38, 'G', sprintf('EPO68\ninoculated');
        35, 'G', sprintf('EPO67\nmock-inoculated');
        29, 'G', sprintf('EPO68\nmock-inoculated')};
dpi_show = [10 12 13 14 16];
[~, j] = ismember(dpi_show, cfg.dpi);
dates_show = cfg.dates(j);
margin = 5;
stretch = 4;   % horizontal stretching factor of the maps (state it in the caption)
nf = size(rows, 1); nd = numel(dpi_show);
% first pass: maps and composites cropped to the leaves (rows: whole
% height; columns: from the first to the last leaf column of the map)
imgs = cell(nf, nd); comps = cell(nf, nd); w = 0;
for f = 1:nf
    for d = 1:nd
        try
            img = rebuild_map_image(rows{f, 1}, dates_show{d}, rows{f, 2}, metadataG, metadataD, ...
                G_pixin1, G_pixin2, G_pixin3, G_pixout, D_pixin1, D_pixin2, D_pixin3, D_pixout, info_table_pix, 1);
            ch = find(any(any(img > 0, 3), 1));
            cols = max(ch(1) - margin, 1):min(ch(end) + margin, size(img, 2));
            imgs{f, d} = img(:, cols, :);
            if strcmpi(rows{f, 2}, 'G')
                comp = composite_image(rows{f, 1}, dates_show{d}, metadataG, ...
                    {G_pixin1, G_pixin2, G_pixin3}, {G_Z1, G_Z2, G_Z3}, band_rgb, size(img, 1), size(img, 2));
            else
                comp = composite_image(rows{f, 1}, dates_show{d}, metadataD, ...
                    {D_pixin1, D_pixin2, D_pixin3}, {D_Z1, D_Z2, D_Z3}, band_rgb, size(img, 1), size(img, 2));
            end
            comps{f, d} = comp(:, cols, :);
            w = max(w, numel(cols));
        catch ME
            fprintf('Pot %d, %s, %d dpi: %s\n', rows{f, 1}, rows{f, 2}, dpi_show(d), ME.message);
        end
    end
end
% one brightness scale for all composites (99th percentile of the leaf
% pixels) and a gamma of 0.6, so that the colours can be compared
allv = cellfun(@(c) c(c > 0), comps, 'UniformOutput', false);
allv = vertcat(allv{:});
sc = prctile(allv, 99);
for k = 1:numel(comps)
    if ~isempty(comps{k})
        comps{k} = min(comps{k} / sc, 1) .^ 0.6;
    end
end
% second pass: all panels padded to the same width (black background)
figure('Color', 'w', 'Position', [50 50 230 * nd 330 * 2 * nf]);
tl = tiledlayout(2 * nf, nd, 'TileSpacing', 'compact', 'Padding', 'compact');
for part = 1:2
    if part == 1, P = imgs; else, P = comps; end
    for f = 1:nf
        for d = 1:nd
            ax = nexttile(tl);
            if isempty(P{f, d})
                text(ax, 0.5, 0.5, 'no image', 'HorizontalAlignment', 'center');
            else
                img = P{f, d};
                pad = w - size(img, 2);
                img = [zeros(size(img, 1), floor(pad / 2), 3), img, zeros(size(img, 1), ceil(pad / 2), 3)];
                % leaves are narrow vertical strips: horizontal stretching
                % (nearest neighbour, no new colours) so that they can be seen
                image(ax, repelem(img, 1, stretch, 1));
                axis(ax, 'image');
            end
            set(ax, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
            if f == 1
                title(ax, sprintf('%d dpi', dpi_show(d)), 'FontSize', 12);
                if d == 1
                    text(ax, -0.03, 1.02, sprintf('(%s)', char('a' + part - 1)), 'Units', 'normalized', ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                        'FontSize', 13, 'FontWeight', 'bold');
                end
            end
            if d == 1
                ylabel(ax, sprintf('%s\n(pot %d)', rows{f, 3}, rows{f, 1}), 'Color', 'k', 'FontSize', 11);
                ax.YLabel.Visible = 'on';
            end
        end
    end
end
exportgraphics(gcf, fullfile(cfg.output_folder, sprintf('Fig6_maps_%s.png', run_mode)), 'Resolution', 300);
close(gcf);
end

function comp = composite_image(pot, date, meta, pixin, Z, band_rgb, nRows, nCols)
% Colour composite (reflectance, not scaled) of the leaf pixels of one image
idx = find(strcmp({meta.fecha}, date) & [meta.planta] == pot, 1);
if isempty(idx)
    error('Pot %d, date %s not found.', pot, date);
end
r = meta(idx).fila; c = meta(idx).columna;
comp = zeros(nRows * nCols, 3);
for i = 1:3
    pix = pixin{i}{r, c};
    if isempty(pix), continue; end
    S = Z{i}{r, c};
    for k = 1:3
        comp(pix, k) = mean(S(:, band_rgb{k}), 2);
    end
end
comp = reshape(max(comp, 0), nRows, nCols, 3);
end
