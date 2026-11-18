function C = project_spectra(D, S, method)
%PROJECT_SPECTRA  MCR contributions of spectra on fixed profiles (Eq. 6)
%
%   C = project_spectra(D, S, method)
%
%   D      : pixels (or observations) x wavelengths
%   S      : MCR profiles (nc x wavelengths)
%   method : 'pinv' -> C = D * pinv(S) (unconstrained least squares)
%            'nnls' -> least squares with C >= 0 (equivalent to fnnls),
%                      consistent with the non-negativity of MCR-ALS

switch method
    case 'pinv'
        C = D * pinv(S);
    case 'nnls'
        C = nnls_gram(S * S', S * D')';
    otherwise
        error('Unknown projection method: %s', method);
end
end
