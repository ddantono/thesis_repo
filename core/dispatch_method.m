% core/dispatch_method.m — updated

function results = dispatch_method(Y, cfg)

    switch lower(cfg.method)

        case 'lasso_enet_var'
            raw = run_lasso_enet_var(Y, cfg);

        case 'mbts_var'
            raw = run_mbts_var(Y, cfg);

        case 'adaptive_lasso_var'
            raw = run_adaptive_lasso_var(Y, cfg);

        case 'ncpen_var'
            raw = run_ncpen_var(Y, cfg);

        case 'hlag_oo_var'
            raw = run_hlag_var(Y, cfg);

        case 'hlag_c_var'
            raw = run_hlag_var(Y, cfg);

        case 'bgr_var'
            raw = run_bgr_var(Y, cfg);

        case 'pds_lm_var'
            raw = run_pds_lm_var(Y, cfg);
            
        case 'msvar'                        
            raw = run_msvar(Y, cfg);

        case 'gglasso_var'                 
            raw = run_gglasso_var(Y, cfg);
 

        otherwise
            error('[dispatch] Unknown method: "%s"', cfg.method);
    end

    results = raw;
end