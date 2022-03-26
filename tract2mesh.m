function [V, F, C] = tract2mesh(varargin)
    
    %% parse inputs
    
    p = inputParser;
    
    % main inputs
    addParameter(p, 'streamlines', []);     % streamlines provided as a cell array - synthetic if not provided
    addParameter(p, 'radius', 0.1);         % individual streamline radius
    addParameter(p, 'vertices', 6);         % number of vertices at cross-section (min: 3)
    addParameter(p, 'colours', []);         % colour-coding: DEC, random, Nx3 matrix (colour per streamline). Returns streamlines indices if empty
    addParameter(p, 'centre', true);        % places middle of tract at the origin
    addParameter(p, 'cores', []);           % number of cores for parallel processing. No parallel processing if empty
    
    % synthetic bundles only - all optional
    addParameter(p, 'nsl', []);             % number of streamlines
    addParameter(p, 'step_size', 1);        % step size    
    
    parse(p, varargin{:});
    fn = fieldnames(p.Results);
    for i = 1:numel(fn)
        eval([fn{i} ' = p.Results.' fn{i} ';']);
    end
    
    % send to parfor branch if desired
    if ~isempty(cores)
        [V, F, C] = tract2mesh_parfor(varargin{:});
        return
    end
    
    %% prepare
    
    r = radius; nv = vertices;
    clear radius vertices
    
    % create synthetic streamlines if needed
    if isempty(streamlines)
        streamlines = simulate_streamline_bundle(nsl, step_size);
    end
    
    % centre streamlines if needed
    if centre
        BB = cellfun(@(x) [max(x) min(x)], streamlines, 'un', 0);
        BBM = cell2mat(BB');
        centroid = (max(BBM(:, 1:3)) + min(BBM(:, 4:6))) / 2;
        streamlines = cellfun(@(x) x - repmat(centroid, [size(x, 1) 1]), streamlines, 'un', 0);
    end
    
    if ~isempty(colours) && isnumeric(colours) && size(colours, 2) == 1
        colours = repmat(colours, [1 3]);
    end
    
    %% convert
    
    nsl = numel(streamlines); % number of streamlines
    [V, F, C] = deal(cell(nsl, 1));       
    
    grad_diff = @(x) norm_vec([gradient(x(:, 1)) gradient(x(:, 2)) gradient(x(:, 3))]);
    
    % generate mesh for individual streamlines
    for i = 1:numel(streamlines)
       
        sl = streamlines{i};
        n_pts = size(sl, 1);                
        
        % gradient vectors - used for rotation and colour-coding
        rtv = grad_diff(sl);
        
        % upwards vectors - avoid twists
        uwv = repmat([0 0 1], [n_pts 1]);
        
        % make uwv perpendicular
        ppv = norm_vec(uwv - repmat(dot(uwv, rtv, 2), [1 3]) .* rtv);
        
        % populate "edges" and "faces" of the streamlines
        for j = 1:nv
        
            % rotation angle
            alpha = 2 * pi * j / nv; 
            
            % Rodrigues formula
            edge_v = ppv .* cos(alpha) + cross(rtv, ppv, 2) .* sin(alpha) + ...
                rtv .* repmat(dot(rtv, ppv, 2), [1 3]) .* (1 - cos(alpha));            
            
            % add streamline vertices + their perpendicular vectors
            % rotated by alpha scaled by radius
            V{i} = [V{i}; sl + edge_v * r]; 
            
            % fill with faces            
            av = (1:n_pts - 1) + n_pts * (j - 1);
            a = [av; av + n_pts + 1; av + n_pts];
            b = [av; av + 1; av + n_pts + 1];
            F{i} = [F{i}; a'; b'];
            
        end
        
        % close the sides
        F{i} = mod(F{i} - 1, size(V{i}, 1)) + 1;
        
        % put the lids on        
        f1 = [ones(1, nv - 2); (1:nv - 2) * n_pts + 1; (2:nv - 1) * n_pts + 1]';
        f2 = f1 + n_pts - 1;
        F{i} = [F{i}; f1; f2(:, [2 1 3])];
        
        % sort out colours
        if nargout > 2 
            if ~isempty(colours) && (ischar(colours) && strcmp(colours, 'DEC') || iscell(colours))
                v_idx = repmat(1:n_pts, [1 nv]);
                v_idx = v_idx(:);
                if iscell(colours)
                    C{i} = colours{i}(v_idx);
                else
                    C{i} = abs(rtv(v_idx, :));
                end
            else
                C{i} = ones(n_pts * nv, 1) * i;
            end
        end
    end
    
    nvc = cellfun(@(x) size(x, 1), V); % number of vertices per cell
    cnvc = [0; cumsum(nvc)]; % cumulative minus the first
    nfc = cellfun(@(x) size(x, 1), F); % number of faces per cell
    cnfc = [0; cumsum(nfc)];
    V = cell2mat(V);    
    F = cell2mat(F);
    C = cell2mat(C);    
    for i = 1:nsl
        F(cnfc(i) + 1:cnfc(i + 1), :) = F(cnfc(i) + 1:cnfc(i + 1), :) + cnvc(i);
    end
    
    % sort out colours
    if nargout > 2        
        if ~isempty(colours) && ischar(colours) && strcmp(colours, 'random')            
            slc = rand(nsl, 3);
            C = slc(C, :);
        elseif ~isempty(colours) && isnumeric(colours) && size(colours, 1) == length(streamlines)
            C = colours(C, :);
        end            
    end
    
    F = F(:, [2 1 3]);
end

function v = norm_vec(v)

    % normalise matrix of vectors
    v = v ./ repmat(sqrt(sum(v .^ 2, 2)), [1 3]);

end