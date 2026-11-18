function [tab, boot_dates] = detection_metrics(R, n_boot, seed, alpha)
%DETECTION_METRICS  Sensitivity, specificity, FPR and balanced accuracy with 95 % CI
%
%   [tab, boot_dates] = detection_metrics(R, n_boot, seed, alpha)
%
%   R : one row per leaf x date observation (out of sample) with columns
%       genotype, inoculated, dpi, pot (pot) and alarm.
%
%   For each genotype and dpi:
%     sensitivity       = proportion of inoculated leaves in alarm
%     specificity       = proportion of mock-inoculated leaves NOT in alarm
%     FPR               = 1 - specificity
%     balanced accuracy = (sensitivity + specificity) / 2
%
%   95 % CI: bootstrap of POTS (the experimental unit), stratified by
%   genotype x treatment: in each resample pots are drawn with replacement
%   within each group and the metrics are recomputed. This accounts for the
%   dependence between leaves of a pot and between dates of a leaf.
%
%   boot_dates : bootstrap distribution of the detection date
%                 (n_boot x number of genotypes).

genos = unique(R.genotype);
dpis  = unique(R.dpi);
nG = numel(genos); nD = numel(dpis);

est = compute(R, genos, dpis);

% --- Pot bootstrap ---
rng(seed);
[grp, ~] = findgroups(R.genotype, R.inoculated);
pots_grp = splitapply(@(p) {unique(p)}, R.pot, grp);
[pots, ~, im] = unique(R.pot);
rows_pot = accumarray(im, (1:height(R))', [], @(x) {x});

bs_sens = nan(n_boot, nG, nD);
bs_spec  = nan(n_boot, nG, nD);
boot_dates = nan(n_boot, nG);
for b = 1:n_boot
    idx = [];
    for g = 1:numel(pots_grp)
        m = pots_grp{g};
        chosen = m(randi(numel(m), numel(m), 1));
        for q = chosen(:)'
            idx = [idx; rows_pot{pots == q}]; %#ok<AGROW>
        end
    end
    Rb = R(idx, :);
    e = compute(Rb, genos, dpis);
    bs_sens(b, :, :) = e.sens;
    bs_spec(b, :, :)  = e.spec;
    for g = 1:nG
        boot_dates(b, g) = detection_date(Rb, genos{g}, alpha);
    end
end
bs_ba = (bs_sens + bs_spec) / 2;

% --- Output table ---
rows = {};
for g = 1:nG
    for d = 1:nD
        ci = @(x) prctile(x(:, g, d), [2.5 97.5]);
        s_ci = ci(bs_sens); e_ci = ci(bs_spec); b_ci = ci(bs_ba);
        rows(end+1, :) = {genos{g}, dpis(d), est.ni(g, d), est.nc(g, d), ...
            est.sens(g, d), s_ci(1), s_ci(2), ...
            est.spec(g, d),  e_ci(1), e_ci(2), ...
            1 - est.spec(g, d), 1 - e_ci(2), 1 - e_ci(1), ...
            est.ba(g, d),   b_ci(1), b_ci(2)}; %#ok<AGROW>
    end
end
tab = cell2table(rows, 'VariableNames', {'genotype','dpi','n_inoc','n_mock', ...
    'sensitivity','sens_CI_low','sens_CI_high', ...
    'specificity','spec_CI_low','spec_CI_high', ...
    'FPR','FPR_CI_low','FPR_CI_high', ...
    'balanced_accuracy','BA_CI_low','BA_CI_high'});
end

% -------------------------------------------------------------------------
function e = compute(R, genos, dpis)
nG = numel(genos); nD = numel(dpis);
e.sens = nan(nG, nD); e.spec = nan(nG, nD);
e.ni = zeros(nG, nD); e.nc = zeros(nG, nD);
for g = 1:nG
    for d = 1:nD
        s  = strcmp(R.genotype, genos{g}) & R.dpi == dpis(d);
        ai = R.alarm(s &  R.inoculated);
        ac = R.alarm(s & ~R.inoculated);
        e.ni(g, d) = numel(ai);
        e.nc(g, d) = numel(ac);
        if ~isempty(ai), e.sens(g, d) = mean(ai); end
        if ~isempty(ac), e.spec(g, d)  = 1 - mean(ac); end
    end
end
e.ba = (e.sens + e.spec) / 2;
end
