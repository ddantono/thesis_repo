function method_cfgs = define_method_configs(pmax, K)
%DEFINE_METHOD_CONFIGS  Return all method configurations for a given pmax.
%
%  method_cfgs = define_method_configs(pmax)
%  method_cfgs = define_method_configs(pmax, K)
%
%  Returns a struct array with 14 entries (one per method/variant).
%  Each entry has:
%    .method   char    framework method name (or 'full_var' for baseline)
%    .label    char    short label for results tables and struct field names
%                      NOTE: labels use underscores only — NO hyphens.
%    .cfg      struct  framework config passed to run_framework
%
%  METHOD LIST:
%    1  Full       Full VAR baseline (no framework call, all-ones support)
%    2  LASSO      LASSO (alpha=1)
%    3  ElasticNet Elastic Net (alpha=0.5)
%    4  mBTS       mBTS — primary reference method from Siggiridou (2016)
%    5  AdapLASSO  Adaptive LASSO
%    6  SCAD       SCAD nonconvex penalty
%    7  MCP        MCP nonconvex penalty
%    8  HLAG_OO    HLAG-OO (dense — benchmark)
%    9  HLAG_C     HLAG-C  (dense — benchmark)
%   10  BGR        BGR/Ridge (dense — benchmark)
%   11  PDS_LM     Post-Double-Selection LM
%   12  msVAR      msVAR (spectral + FDR)
%   13  gLASSO_V   Group LASSO (grouped by variable)
%   14  gLASSO_L   Group LASSO (grouped by lag)

if nargin < 2 || isempty(K)
    K = [];
end

method_cfgs = struct('method', {}, 'label', {}, 'cfg', {});
n = 0;

% 1. Full VAR baseline
n = n+1;
method_cfgs(n).method = 'full_var';
method_cfgs(n).label  = 'Full';
method_cfgs(n).cfg    = struct('method','full_var','pmax',pmax);

% 2. LASSO
n = n+1;
method_cfgs(n).method = 'lasso_enet_var';
method_cfgs(n).label  = 'LASSO';
method_cfgs(n).cfg    = struct('method','lasso_enet_var','pmax',pmax,'alpha',1.0);

% 3. Elastic Net
n = n+1;
method_cfgs(n).method = 'lasso_enet_var';
method_cfgs(n).label  = 'ElasticNet';
method_cfgs(n).cfg    = struct('method','lasso_enet_var','pmax',pmax,'alpha',0.5);

% 4. mBTS  (primary reference method from Siggiridou 2016)
n = n+1;
method_cfgs(n).method = 'mbts_var';
method_cfgs(n).label  = 'mBTS';
method_cfgs(n).cfg    = struct('method','mbts_var','pmax',pmax);

% 5. Adaptive LASSO
n = n+1;
method_cfgs(n).method = 'adaptive_lasso_var';
method_cfgs(n).label  = 'AdapLASSO';
method_cfgs(n).cfg    = struct('method','adaptive_lasso_var','pmax',pmax);

% 6. SCAD
n = n+1;
method_cfgs(n).method = 'ncpen_var';
method_cfgs(n).label  = 'SCAD';
method_cfgs(n).cfg    = struct('method','ncpen_var','pmax',pmax,'penalty','SCAD');

% 7. MCP
n = n+1;
method_cfgs(n).method = 'ncpen_var';
method_cfgs(n).label  = 'MCP';
method_cfgs(n).cfg    = struct('method','ncpen_var','pmax',pmax,'penalty','MCP');

% 8. HLAG-OO  (dense output — OLS-CGCI on full support = full-VAR CGCI)
n = n+1;
method_cfgs(n).method = 'hlag_oo_var';
method_cfgs(n).label  = 'HLAG_OO';
method_cfgs(n).cfg = struct('method','hlag_oo_var','pmax',pmax,'struct','HLAGOO');

% 9. HLAG-C  (dense output)
n = n+1;
method_cfgs(n).method = 'hlag_c_var';
method_cfgs(n).label  = 'HLAG_C';
method_cfgs(n).cfg = struct('method','hlag_c_var','pmax',pmax,'struct','HLAGC');

% 10. BGR/Ridge  (dense output)
n = n+1;
method_cfgs(n).method = 'bgr_var';
method_cfgs(n).label  = 'BGR';
method_cfgs(n).cfg    = struct('method','bgr_var','pmax',pmax);

% 11. PDS-LM
n = n+1;
method_cfgs(n).method = 'pds_lm_var';
method_cfgs(n).label  = 'PDS_LM';
method_cfgs(n).cfg    = struct('method','pds_lm_var','pmax',pmax);

% 12. msVAR
n = n+1;
method_cfgs(n).method = 'msvar';
method_cfgs(n).label  = 'msVAR';
method_cfgs(n).cfg    = struct('method','msvar','p_seq',1:pmax);

% 13. Group LASSO — variable grouping
n = n+1;
method_cfgs(n).method = 'gglasso_var';
method_cfgs(n).label  = 'gLASSO_V';
method_cfgs(n).cfg    = struct('method','gglasso_var','pmax',pmax,'group_by','variable');

% 14. Group LASSO — lag grouping
n = n+1;
method_cfgs(n).method = 'gglasso_var';
method_cfgs(n).label  = 'gLASSO_L';
method_cfgs(n).cfg    = struct('method','gglasso_var','pmax',pmax,'group_by','lag');

end
