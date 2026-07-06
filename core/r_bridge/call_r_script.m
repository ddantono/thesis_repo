% core/r_bridge/call_r_script.m

function [success, log_msg] = call_r_script(r_exe, r_script, varargin)
% CALL_R_SCRIPT  Executes any R script via system call.
%
%   [success, log_msg] = call_r_script(r_exe, r_script, arg1, arg2, ...)
%
%   All varargin are passed as quoted command-line arguments to R.
%   R script receives them via commandArgs(trailingOnly=TRUE).
%
%   The last varargin is assumed to be the expected output file path
%   and is checked for existence after the call.
%
%   USAGE
%     call_r_script('Rscript', 'method.R', in_file, cfg_file, out_file)

    % Build argument string
    arg_str = '';
    for i = 1:numel(varargin)
        arg_str = [arg_str, sprintf(' "%s"', varargin{i})]; %#ok<AGROW>
    end

    cmd = sprintf('"%s" "%s"%s', r_exe, r_script, arg_str);

    [status, log_msg] = system(cmd);
    log_msg           = strtrim(log_msg);

    % Expected output file is the last argument
    expected_out = varargin{end};
    file_exists  = exist(expected_out, 'file') == 2;

    if status ~= 0
        success = false;
        if isempty(log_msg)
            log_msg = sprintf('R exited with status %d.', status);
        end
    elseif ~file_exists
        success = false;
        log_msg = sprintf('R completed (status=0) but output file not found: %s', ...
                          expected_out);
    else
        success = true;
        log_msg = '';
    end
end