function X = nnls_gram(G, B)
%NNLS_GRAM  Non-negative least squares for many columns at once
%
%   X = nnls_gram(G, B)
%
%   For each column j, solves  min ||y_j - A x_j||^2  subject to x_j >= 0,
%   from G = A'*A (k x k) and B = A'*Y (k x N). Returns X (k x N).
%
%   For small k (4-5 MCR components) it is EXACT: it tries the 2^k - 1
%   possible supports, computes the least-squares solution on each, and
%   keeps for each column the feasible (all >= 0) solution with the lowest
%   error. Gives the same result as fnnls, vectorised over millions of pixels.
%
%   Objective (without the constant y'y):  f(x) = x'Gx - 2 b'x

k = size(G, 1);
N = size(B, 2);
X = zeros(k, N);
best_f = zeros(1, N);          % x = 0 is always feasible, with f = 0

for s = 1:(2^k - 1)
    sop = logical(bitget(s, 1:k));
    Xs = zeros(k, N);
    Xs(sop, :) = G(sop, sop) \ B(sop, :);
    feasible = all(Xs(sop, :) >= 0, 1);
    if ~any(feasible), continue; end
    f = sum(Xs .* (G * Xs), 1) - 2 * sum(B .* Xs, 1);
    improves = feasible & (f < best_f);
    X(:, improves) = Xs(:, improves);
    best_f(improves) = f(improves);
end
end
