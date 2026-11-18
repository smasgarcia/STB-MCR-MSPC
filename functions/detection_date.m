function [dpi_det, tabp] = detection_date(R, genotype, alpha)
%DETECTION_DATE  First dpi with statistically significant leaf-level detection
%
%   [dpi_det, tabp] = detection_date(R, genotype, alpha)
%
%   Criterion: first dpi >= 0 at which the proportion of
%   inoculated leaves in alarm is higher than that of mock-inoculated
%   leaves of the same genotype (one-sided Fisher exact test, Holm
%   correction over all evaluated dpi).
%
%   R       : table with genotype, inoculated, dpi and alarm (logical)
%   dpi_det : first significant dpi (NaN if none)
%   tabp    : table with p and Holm-adjusted p per dpi
%
%   The test treats leaves as independent within each dpi; the uncertainty
%   of the date is therefore also given by a pot bootstrap (see s03).

sel  = strcmp(R.genotype, genotype) & R.dpi >= 0;
dpis = unique(R.dpi(sel));
p    = nan(numel(dpis), 1);
prop_inoc = nan(numel(dpis), 1);
prop_mock = nan(numel(dpis), 1);

for d = 1:numel(dpis)
    s  = sel & R.dpi == dpis(d);
    ai = R.alarm(s &  R.inoculated);
    ac = R.alarm(s & ~R.inoculated);
    if isempty(ai) || isempty(ac), continue; end
    prop_inoc(d) = mean(ai);
    prop_mock(d) = mean(ac);
    tabla = [sum(ai) sum(~ai); sum(ac) sum(~ac)];
    [~, p(d)] = fishertest(tabla, 'Tail', 'right');
end

padj = holm(p);
k = find(padj < alpha, 1);
if isempty(k)
    dpi_det = NaN;
else
    dpi_det = dpis(k);
end
tabp = table(dpis, prop_inoc, prop_mock, p, padj, ...
    'VariableNames', {'dpi','prop_alarm_inoc','prop_alarm_mock','p','p_holm'});
end
