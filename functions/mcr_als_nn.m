function [C, S, info] = mcr_als_nn(D, S0, opts)
%MCR_ALS_NN  MCR-ALS with non-negativity (C and S) and unit-length spectra
%
%   [C, S, info] = mcr_als_nn(D, S0, opts)
%
%   D    : pixels x wavelengths
%   S0   : initial spectra (nc x wavelengths), e.g. from SIMPLISMA
%   opts  : .max_iter  maximum number of iterations
%          .tol       convergence: relative change in lack of fit (%)
%
%   Constraints: non-negativity of concentrations and spectra (exact NNLS,
%   see nnls_gram) and spectra normalised to unit length. Implemented here
%   so that the model can be refitted within the cross-validation.
%
%   info.lof         lack of fit (%) = 100*sqrt(SSE/SST)
%   info.r2          explained variance (%) = 100*(1 - SSE/SST)
%   info.n_iter      number of iterations
%   info.converged  true if the tolerance was reached

S = S0 ./ vecnorm(S0, 2, 2);
sst = sum(D(:).^2);
lof_prev = Inf;
info.converged = false;

for it = 1:opts.max_iter
    % C step: for each pixel, min ||d - c S|| with c >= 0
    C = nnls_gram(S * S', S * D')';
    % S step: for each wavelength, min ||D(:,j) - C s_j|| with s_j >= 0
    S = nnls_gram(C' * C, C' * D);
    % Spectra normalised to unit length (scale moved to C)
    nrm = vecnorm(S, 2, 2);
    nrm(nrm == 0) = 1;
    S = S ./ nrm;
    C = C .* nrm';
    % Error without building the full residual matrix (saves memory)
    sse = sst - 2 * sum(sum(C .* (D * S'))) + sum(sum((C' * C) .* (S * S')));
    lof = 100 * sqrt(max(sse, 0) / sst);
    if abs(lof_prev - lof) / lof * 100 < opts.tol
        info.converged = true;
        break;
    end
    lof_prev = lof;
end

info.lof = lof;
info.r2 = 100 * (1 - max(sse, 0) / sst);
info.n_iter = it;
end
