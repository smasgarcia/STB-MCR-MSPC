function run_mode = mode_name(cfg)
%MODE_NAME  Suffix of the output files according to the config options
%
%   'full'   default options (MCR-ALS on all calibration leaves)
%   'oneclass'   MCR-ALS on mock-inoculated calibration leaves only
%   ..._initfull  MCR-ALS initialised with the full-model profiles
%   ..._refgeno  genotype-specific MSPC reference
%   ..._nnls     non-negative projection
%   ..._k4       number of PCs set by the user (here 4)

if cfg.mcr_mock_only, run_mode = 'oneclass'; else, run_mode = 'full'; end
if isfield(cfg, 'mcr_init') && strcmp(cfg.mcr_init, 'full_model')
    run_mode = [run_mode '_initfull'];
end
if isfield(cfg.mspc, 'ref_by_genotype') && cfg.mspc.ref_by_genotype
    run_mode = [run_mode '_refgeno'];
end
if isfield(cfg, 'projection') && strcmp(cfg.projection, 'nnls')
    run_mode = [run_mode '_nnls'];
end
if isfield(cfg.mspc, 'k') && ~isempty(cfg.mspc.k)
    run_mode = sprintf('%s_k%d', run_mode, cfg.mspc.k);
end
end
