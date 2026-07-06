function cfg = mbts_var_config()
% Default config for mBTS-VAR wrapper.

    cfg.method   = 'mbts_var';
    cfg.p        = [];    % [] → use cfg.pmax
    cfg.pmax     = 2;    % maximum lag order to search

    % --- output ----------------------------------------------------------
    cfg.save     = true;
    cfg.verbose  = false;
end