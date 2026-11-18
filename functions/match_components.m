function [S_ord, ord_idx, r] = match_components(S, Sref)
%MATCH_COMPONENTS  Reorder components to match a reference
%
%   [S_ord, ord_idx, r] = match_components(S, Sref)
%
%   MCR-ALS does not guarantee the order of the components, so each fit
%   (each cross-validation fold) may return them permuted. The permutation
%   that maximises the total correlation with the reference profiles
%   (full model, soptTots) is chosen. Used only to LABEL the components; it
%   plays no role in the fit and introduces no information leakage.
%
%   r : correlation of each reordered component with its reference
%       (measure of profile stability).

nc = size(S, 1);
assert(size(Sref, 1) == nc, 'S and Sref must have the same number of components.');

R = corr(S', Sref');                      % nc x nc
P = perms(1:nc);
points = zeros(size(P, 1), 1);
for i = 1:size(P, 1)
    points(i) = sum(R(sub2ind([nc nc], P(i, :), 1:nc)));
end
[~, ib] = max(points);
ord_idx = P(ib, :);
S_ord = S(ord_idx, :);
r = diag(corr(S_ord', Sref'))';
end
