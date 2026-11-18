function M = mspc_fit(X, opts)
%MSPC_FIT  PCA reference model and T2 / Q control limits
%
%   M = mspc_fit(X, opts)
%
%   X    : reference observations (mock-inoculated CALIBRATION leaves) x
%          variables (mean MCR-ALS contributions)
%   opts  : .alpha    significance level (0.05)
%          .var_pca  fraction of variance used to choose k (if .k is empty)
%          .k        number of PCs set by the user (optional)
%
%   Uses svd instead of pca() to avoid conflicts with other toolboxes.
%
%   T2 limits computed (one is chosen in config):
%     T2lim_chi2 : chi2inv(1-alpha, k)
%     T2lim_F1   : k(N-1)/(N-k) F(k, N-k)       (phase I)
%     T2lim_F    : k(N^2-1)/(N(N-k)) F(k, N-k)  (phase II: NEW observations,
%                  appropriate for validation leaves)
%   Q (SPE) limit: Jackson and Mudholkar (1979), from the eigenvalues of the
%   PCs NOT retained. If k = number of variables, Q is always 0 and the
%   limit is NaN.

M.mean_vec = mean(X, 1);
Xc = X - M.mean_vec;
N  = size(X, 1);
[~, Sv, V] = svd(Xc, 'econ');
lambda = diag(Sv).^2 / (N - 1);
M.expl = 100 * lambda / sum(lambda);

if isfield(opts, 'k') && ~isempty(opts.k)
    k = opts.k;
else
    k = find(cumsum(M.expl) >= 100 * opts.var_pca - 1e-9, 1);
end

M.k = k;
M.N = N;
M.P = V(:, 1:k);
M.lambda = lambda(1:k);
M.lambda_all = lambda;
alpha = opts.alpha;

% --- T2 limits ---
M.T2lim_chi2 = chi2inv(1 - alpha, k);
M.T2lim_F1   = k * (N - 1) / (N - k) * finv(1 - alpha, k, N - k);
M.T2lim_F    = k * (N^2 - 1) / (N * (N - k)) * finv(1 - alpha, k, N - k);

% --- Q limit (Jackson-Mudholkar) ---
lr = lambda(k+1:end);
if isempty(lr) || sum(lr) <= eps
    M.Qlim = NaN;
else
    t1 = sum(lr); t2 = sum(lr.^2); t3 = sum(lr.^3);
    h0 = 1 - 2 * t1 * t3 / (3 * t2^2);
    ca = norminv(1 - alpha);
    M.Qlim = t1 * (ca * sqrt(2 * t2 * h0^2) / t1 + 1 + t2 * h0 * (h0 - 1) / t1^2)^(1 / h0);
end

% --- Empirical limit (percentile of the reference itself) ---
[T2ref, Qref] = mspc_apply(M, X);
M.T2lim_emp = prctile(T2ref, 100 * (1 - alpha));
M.Qlim_emp  = prctile(Qref, 100 * (1 - alpha));
end
