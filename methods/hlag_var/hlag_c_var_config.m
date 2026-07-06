% methods/hlag_var/hlag_c_var_config.m

function cfg = hlag_c_var_config()
% Convenience config for HLAG Componentwise VAR.

    cfg        = hlag_var_config();
    cfg.method = 'hlag_c_var';
    cfg.struct = 'HLAGC';
end