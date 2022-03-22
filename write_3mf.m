function write_3mf(varargin)

    %% parse inputs
    
    p = inputParser;
    
    % main inputs
    addRequired(p, 'fname');            % output file name
    addRequired(p, 'V');                % vertices
    addRequired(p, 'F');                % faces
    addParameter(p, 'C', []);           % colours
    addParameter(p, 'N', false);        % normals
    addParameter(p, 'precision', 3);    % precision digits (vertex coordinates)
    
    parse(p, varargin{:});
    fn = fieldnames(p.Results);
    for i = 1:numel(fn)
        eval([fn{i} ' = p.Results.' fn{i} ';']);
    end
    
    %% prepare
    
    % convert to HEX
    if ~isempty(C)
        H(:, 1:6) = reshape(sprintf('%02X', round(C * 255)'), 6, [])'; 
        C = H;
    end
    
    if N
        N = calculate_normals(V, F);
    end
    
    nv = num2str(size(V, 1));
    nf = num2str(size(F, 1));
    pr = '32';
    
    
    
    
    
    
    % Temporary files
    files = tempFiles(filename); 
    writeTempFiles(files);   
    
    % 3D model
    write3DModel(files , vertices , faces, colors);

    % Builds 3mf file
    package3mf(files)

    % Remove temp files
    cleanTempFiles(files)   
    
end

function N = calculate_normals(V, F)

    nf = size(F, 1);
    nv = size(V, 1);
    N = zeros(nv, 3);
    d = @(v) sqrt(sum(v .^ 2, 2));    
    
    % face normals
    NF = cross(V(F(:, 2), :) - V(F(:, 1), :), V(F(:, 3), :) - V(F(:, 1), :)); 
    D = d(NF); 
    D (D < eps) = 1;
    NF = NF ./ repmat(D, 1, 3);

    % vertex normals
    for i = 1:nf
        f = F(i, :);
        for j = 1:3
            N(f(j), :) = N(f(j), :) + NF(i, :);
        end
    end
    D = d(N);
    D(D < eps) = 1;
    N = N ./ repmat(D, 1, 3);    
end