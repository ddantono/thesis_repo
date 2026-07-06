% methods/hlag_var/hlag_oo_var_config.m

function cfg = hlag_oo_var_config()
% Convenience config for HLAG Own/Other VAR.

    cfg        = hlag_var_config();
    cfg.method = 'hlag_oo_var';
    cfg.struct = 'HLAGOO';
end