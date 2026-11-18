function padj = holm(p)
%HOLM  Holm-Bonferroni adjustment (NaN values are ignored)
%
%   padj = holm(p)

padj = nan(size(p));
ok = ~isnan(p);
q = p(ok);
m = numel(q);
if m == 0, return; end
[qs, i] = sort(q(:));
adj = min(1, cummax((m:-1:1)' .* qs));
tmp = nan(m, 1);
tmp(i) = adj;
padj(ok) = tmp;
end
