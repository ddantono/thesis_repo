function cfg = elastic_net_var_config()
% Convenience config for Elastic Net — thin wrapper over lasso_var_config.
% Only alpha differs from pure LASSO.

    cfg         = lasso_var_config();
    cfg.method  = 'lasso_enet_var';   % same wrapper, same dispatch key
    cfg.alpha   = 0.5;           % equal L1/L2 mix — tune as needed
end