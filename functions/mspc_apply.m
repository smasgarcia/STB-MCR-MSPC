function [T2, Q, Tsc] = mspc_apply(M, X)
%MSPC_APPLY  Project observations onto the PCA reference model
%
%   [T2, Q, Tsc] = mspc_apply(M, X)
%
%   X   : observations x variables (same variables as in mspc_fit)
%   T2  : Hotelling T2 = sum(t.^2 ./ lambda)
%   Q   : squared prediction error (SPE) = ||x - t P'||^2
%   Tsc : PCA scores

Xc  = X - M.mean_vec;
Tsc = Xc * M.P;
T2  = sum(Tsc.^2 ./ M.lambda', 2);
E   = Xc - Tsc * M.P';
Q   = sum(E.^2, 2);
end
