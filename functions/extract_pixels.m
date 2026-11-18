function [D, id_obs] = extract_pixels(result, T, rows, idx_lambda, max_pix)
%EXTRACT_PIXELS  Stack the pixel spectra of several leaf x date observations
%
%   [D, id_obs] = extract_pixels(result, T, rows, idx_lambda, max_pix)
%
%   rows    : rows of T to extract
%   max_pix  : maximum number of pixels per observation (random subsample
%              without replacement). Inf = all pixels.
%   D        : pixels x wavelengths
%   id_obs   : for each pixel, the row of T it comes from
%
%   The subsample uses the global random generator: set rng before calling
%   this function for reproducibility.

rows = rows(:)';
n_tot = sum(min(T.npix(rows), max_pix));
D = zeros(n_tot, numel(idx_lambda));
id_obs = zeros(n_tot, 1);

p = 0;
for r = rows
    block = result{T.i_res(r)}.data(T.row_start(r):T.row_end(r), idx_lambda);
    if size(block, 1) > max_pix
        block = block(randperm(size(block, 1), max_pix), :);
    end
    n = size(block, 1);
    D(p+1:p+n, :) = block;
    id_obs(p+1:p+n) = r;
    p = p + n;
end
end
