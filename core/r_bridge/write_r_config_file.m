% core/r_bridge/write_r_config_file.m

function write_r_config_file(cfg_struct, filepath)
% WRITE_R_CONFIG_FILE  Writes any scalar/string struct to CSV for R.
%
%   Format: two columns — param,value
%   Supports: scalar numeric, logical, char fields.
%   Skips:    struct, cell, matrix fields (not R-serializable here).
%
%   USAGE
%     cfg_r.p           = 2;
%     cfg_r.n_folds     = 10;
%     cfg_r.cv_criterion= 'lambda.min';
%     write_r_config_file(cfg_r, filepath);

    fid = fopen(filepath, 'w');
    if fid == -1
        error('[write_r_config_file] Cannot open: %s', filepath);
    end

    fprintf(fid, 'param,value\n');

    fields = fieldnames(cfg_struct);
    for i = 1:numel(fields)
        f = fields{i};
        v = cfg_struct.(f);

        if ischar(v) || isstring(v)
            fprintf(fid, '%s,%s\n', f, char(v));

        elseif islogical(v) && isscalar(v)
            fprintf(fid, '%s,%d\n', f, int32(v));

        elseif isnumeric(v) && isscalar(v)
            if v == floor(v)
                fprintf(fid, '%s,%d\n', f, int64(v));
            else
                fprintf(fid, '%s,%.10g\n', f, v);
            end

        else
            % Skip non-scalar / struct / cell fields silently
        end
    end

    fclose(fid);
end