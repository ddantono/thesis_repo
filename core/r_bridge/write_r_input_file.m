% core/r_bridge/write_r_input_file.m

function write_r_input_file(Y, filepath)
% WRITE_R_INPUT_FILE  Writes Y matrix to CSV for any R method.
%
%   Format: header y1..yK, then numeric data at full precision.

    [~, K] = size(Y);
    headers     = arrayfun(@(k) sprintf('y%d',k), 1:K, 'UniformOutput', false);
    header_line = strjoin(headers, ',');

    fid = fopen(filepath, 'w');
    if fid == -1
        error('[write_r_input_file] Cannot open: %s', filepath);
    end
    fprintf(fid, '%s\n', header_line);
    fclose(fid);

    dlmwrite(filepath, Y, '-append', 'delimiter', ',', 'precision', 15);
end