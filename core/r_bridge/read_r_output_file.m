% core/r_bridge/read_r_output_file.m

function result = read_r_output_file(filepath)
% READ_R_OUTPUT_FILE  Reads any R output CSV into a MATLAB struct.
%
%   Expected file format — sections separated by named headers:
%
%     <SECTION_NAME>
%     <rows of comma-separated values>
%     <SECTION_NAME>
%     ...
%
%   Numeric sections  → parsed as double matrix
%   Key-value sections→ parsed as struct fields
%
%   Convention used by all R methods in this framework:
%     - Section "B"    : [K*p x K] coefficient matrix
%     - Section "meta" : key,value pairs (scalars and vectors)
%
%   OUTPUT
%     result : struct with one field per section name
%              result.B    — numeric matrix
%              result.meta — struct of key-value pairs

    if ~exist(filepath, 'file')
        error('[read_r_output_file] File not found: %s', filepath);
    end

    % Read all lines
    fid   = fopen(filepath, 'r');
    lines = {};
    tline = fgetl(fid);
    while ischar(tline)
        lines{end+1} = strtrim(tline); %#ok<AGROW>
        tline = fgetl(fid);
    end
    fclose(fid);

    lines = lines(~cellfun(@isempty, lines));

    % Locate section headers — a header is a line with no commas
    % that is immediately followed by data lines
    is_header = false(numel(lines), 1);
    for i = 1:numel(lines)
        if isempty(strfind(lines{i}, ','))
            is_header(i) = true;
        end
    end

    header_idx   = find(is_header);
    header_names = lines(header_idx);

    result = struct();

    for h = 1:numel(header_idx)
        name     = header_names{h};
        i_start  = header_idx(h) + 1;
        if h < numel(header_idx)
            i_end = header_idx(h+1) - 1;
        else
            i_end = numel(lines);
        end

        section_lines = lines(i_start:i_end);
        section_lines = section_lines(~cellfun(@isempty, section_lines));

        % Determine if section is key-value or pure numeric matrix
        first = section_lines{1};
        parts = strsplit(first, ',');
        first_val = str2double(parts{1});

        if isnan(first_val)
            % Key-value section
            result.(name) = i_parse_kv(section_lines);
        else
            % Numeric matrix section
            result.(name) = i_parse_matrix(section_lines);
        end
    end
end


% --- private helpers ---

function s = i_parse_kv(lines)
    s = struct();
    for i = 1:numel(lines)
        parts = strsplit(lines{i}, ',');
        key   = matlab.lang.makeValidName(strtrim(parts{1}));
        vals  = strtrim(parts(2:end));

        nums = str2double(vals);

        if numel(vals) == 1 && isnan(nums(1))
            % Check if literal 'NaN' — store as numeric NaN, not char
            if strcmpi(strtrim(vals{1}), 'NaN')
                s.(key) = NaN;
            else
                % Single string value
                s.(key) = char(vals{1});
            end
        elseif all(isfinite(nums))
            % Numeric — scalar or vector
            s.(key) = nums;
        else
            % Mixed or unparseable
            s.(key) = strjoin(vals, ',');
        end
    end
end
function M = i_parse_matrix(lines)
% Parses lines of comma-separated numbers into a numeric matrix.

    n_rows = numel(lines);
    first  = str2double(strsplit(lines{1}, ','));
    n_cols = numel(first);
    M      = zeros(n_rows, n_cols);
    M(1,:) = first;
    for i = 2:n_rows
        M(i,:) = str2double(strsplit(lines{i}, ','));
    end
end