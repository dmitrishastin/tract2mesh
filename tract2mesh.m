function [V, F, C] = tract2mesh(varargin)
    
    %% parse inputs
    
    p = inputParser;
    
    % main inputs
    addParameter(p, 'streamlines', []);     % streamlines provided as a cell array - synthetic if not provided
    addParameter(p, 'radius', 0.1);         % individual streamline radius
    addParameter(p, 'vertices', 6);         % number of vertices at cross-section (min: 3)
    addParameter(p, 'colours', []);         % colour-coding: DEC, random, Nx3 matrix (colour per streamline). Returns streamlines indices if empty
    addParameter(p, 'centre', true);        % places middle of tract at the origin
    
    % synthetic bundles only - all optional
    addParameter(p, 'nsl', []);             % number of streamlines
    addParameter(p, 'step_size', 1);        % step size    
    
    parse(p, varargin{:});
    fn = fieldnames(p.Results);
    for i = 1:numel(fn)
        eval([fn{i} ' = p.Results.' fn{i} ';']);
    end
    
    %% prepare
    
    r = radius;
    nv = vertices;
    clear radius vertices
    
    % create synthetic streamlines if needed
    if isempty(streamlines)
        streamlines = simulate_streamline_bundle(nsl, step_size);
    end
    
    % centre streamlines if needed
    if centre
        BB = cellfun(@(x) [max(x) min(x)], streamlines, 'un', 0);
        BBM = cell2mat(BB);
        centroid = (max(BBM(:, 1:3)) + min(BBM(:, 4:6))) / 2;
        streamlines = cellfun(@(x) x - repmat(centroid, [size(x, 1) 1]), streamlines, 'un', 0);
    end
    
    if ~isempty(colours) && isnumeric(colours) && size(colours, 2) == 1
        colours = repmat(colours, [1 3]);
    end
    
    %% convert
    
    [V, F, C] = deal([]);    
    anv = 0; % intercept for vertex enumeration (faces array)
    
    % generate cross-sectional disc    
    [dv, df] = gen_disc(nv, r);  
    
    % generate and grow segments
    for i = 1:numel(streamlines)
       
        sl = streamlines{i};
        n_pts = size(sl, 1);
        
        % gradient vectors - used for rotation and colour-coding
        rtv = norm_vec([gradient(sl(:, 1)) gradient(sl(:, 2)) gradient(sl(:, 3))]);
        norot = all(abs(rtv) == repmat([0 0 1], [n_pts 1]), 2);
        
        % pre-generate cross and dot products for rotation matrices
        cross_v = cross(repmat([0 0 1], [n_pts, 1]), rtv, 2);
        dot_v = dot(repmat([0 0 1], [n_pts, 1]), rtv, 2);
        
        % position the first disc
        if ~norot(1)                                    
            R = rot_mtx(cross_v(1, :), dot_v(1, :), sl(1, :));
            v1 = (R * [dv ones(nv, 1)]')';            
        else
            v1 = dv + sl(1, :);
        end
        V = [V; v1(:, 1:3)];
        F = [F; df + anv];
       
        for j = 2:n_pts
            
            % position the vertices of the next disc             
            if ~norot(j)
                R = rot_mtx(cross_v(j, :), dot_v(j, :), sl(j, :));
                v1 = (R * [dv ones(nv, 1)]')';
            else
                v1 = dv + sl(j, :);
            end 
            
            V = [V; v1(:, 1:3)];
            
            % generate the walls
            av = 1:nv;
            a = [av; av+1; av+nv+1]';
            b = [av; av+nv+1; av+nv]';
            a(end, :) = [nv 1 nv+1];
            b(end, :) = [nv nv+1 nv+nv];
            
            F = [F; a + anv; b + anv];
            anv = anv + nv;
            
        end
        
        % final disc        
        F = [F; df + anv];
        anv = anv + nv;
        
        % sort out colours
        if nargout > 2           
            if ~isempty(colours) && ischar(colours)
                switch colours
                    case 'DEC'
                        v_idx = repmat(1:n_pts, [nv 1]);
                        v_idx = v_idx(:);
                        C = [C; abs(rtv(v_idx, :))];
                    case 'random'
                        slc = rand(1, 3);
                        C = [C; repmat(slc, [n_pts * nv 1])];
                end
            elseif ~isempty(colours) && isnumeric(colours) && size(colours, 1) == length(streamlines)
                C = [C; repmat(colours(i, :), [n_pts * nv 1])];
            else                 
                C = [C; ones(n_pts * nv, 1) * i];
            end            
        end
    end    
end

function [v, f] = gen_disc(nf, r)

    % generate a mesh representation of a disc
    alpha = 2 * pi / nf; 
    v = [sin(alpha .* (0:nf - 1)); cos(alpha .* (0:nf - 1)); zeros(1, nf)]';
    v = (v - mean(v)) * r;
    f = [ones(1, nf - 2); 2:nf - 1; 3:nf]';        

end

function R = rot_mtx(cv, dv, O)

    % generate transform matrix for each disc
    ssc = [0 -cv(3) cv(2); cv(3) 0 -cv(1); -cv(2) cv(1) 0];
    R = eye(3) + ssc + ssc ^ 2 * (1 - dv) / sqrt(sum(cv .^ 2)) ^ 2;
    R = [R(1, :) O(1); R(2, :) O(2); R(3, :) O(3); 0 0 0 1];    

end

function v = norm_vec(v)

    % normalise matrix of vectors
    v = v ./ repmat(sqrt(sum(v .^ 2, 2)), [1 3]);

end