function [S0, sel] = simplisma_spectra(D, nc, offset_pct)
%SIMPLISMA_SPECTRA  Initial spectral estimates (purest pixels)
%
%   [S0, sel] = simplisma_spectra(D, nc, offset_pct)
%
%   D          : pixels x wavelengths (use a subsample of ~20,000 pixels)
%   nc         : number of components
%   offset_pct : SIMPLISMA offset as % of the maximum mean (the "noise
%                level" of the MCR-ALS GUI, typically 1-10 %)
%
%   Applies SIMPLISMA (Windig and Guilment, 1991) with pixels as
%   "variables": selects the nc purest pixels and returns their spectra
%   (S0, nc x wavelengths) and indices (sel). Weights are computed with the
%   Schur complement, without building the full pixel x pixel correlation
%   matrix.

X  = D';                                  % wavelengths x pixels
m  = size(X, 1);
mu = mean(X, 1);
sd = std(X, 1, 1);
f  = offset_pct / 100 * max(mu);

base_purity = sd ./ (mu + f);             % purity of each pixel
Y   = X ./ sqrt(sd.^2 + (mu + f).^2);     % normalised data (Windig)
Rjj = sum(Y.^2, 1) / m;                   % diagonal of the correlation matrix

sel = zeros(1, nc);
for k = 1:nc
    if k == 1
        w = Rjj;
    else
        K   = sel(1:k-1);
        RjK = (Y' * Y(:, K)) / m;         % pixels x (k-1)
        RKK = RjK(K, :);
        % det([Rjj RjK; RKj RKK]) / det(RKK) = Rjj - RjK RKK^-1 RKj
        w = Rjj - sum((RjK / RKK) .* RjK, 2)';
    end
    [~, sel(k)] = max(w .* base_purity);
end

S0 = D(sel, :);
end
