function [V, F, C] = tract2mesh_parfor(varargin)
    
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
    
    %% prepare
    
    r = radius;
    nv = vertices;
    clear radius vertices
        
    PP = gcp('nocreate');
    if isempty(PP) || PP.NumWorkers ~= cores
        delete(PP);
        cust_clust = parcluster; 
        cust_clust.NumWorkers = cores;
        parpool(cust_clust, 'IdleTimeout', 300);
    end
    clear PP    
    
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
    
    % generate cross-sectional disc    
    [dv, df] = gen_disc(nv, r);  
    
    % generate and grow segments
    parfor i = 1:nsl
       
        sl = streamlines{i};
        n_pts = size(sl, 1);
        anv = 0; % intercept for vertex enumeration (faces array)        
        
        % gradient vectors - used for rotation and colour-coding
        rtv = norm_vec([gradient(sl(:, 1)) gradient(sl(:, 2)) gradient(sl(:, 3))]);
        norot = all(abs(rtv) == repmat([0 0 1], [n_pts 1]), 2);
        
        % pre-generate cross and dot products for rotation matrices
        cross_v = cross(repmat([0 0 1], [n_pts, 1]), rtv, 2);
        dot_v = dot(repmat([0 0 1], [n_pts, 1]), rtv, 2);
       
        for j = 1:n_pts
            
            % place disc vertices around each streamline point
            if ~norot(j)
                R = rot_mtx(cross_v(j, :), dot_v(j, :), sl(j, :));
                v1 = (R * [dv ones(nv, 1)]')';
            else
                v1 = dv + sl(j, :);
            end             
            V{i} = [V{i}; v1(:, 1:3)];
            
            % first point: cap off the streamline, no side walls
            if j == 1
                F{i} = df;
                continue
            end
            
            % subseqent points: generate walls, leave cross-section empty
            av = 1:nv;
            a = [av; av+1; av+nv+1]';
            b = [av; av+nv+1; av+nv]';
            a(end, :) = [nv 1 nv+1];
            b(end, :) = [nv nv+1 nv+nv];            
            F{i} = [F{i}; a + anv; b + anv];
            anv = anv + nv;
            
        end
        
        % final point: cap off the streamline
        F{i} = [F{i}; df + anv];
        
        % sort out colours
        if ~isempty(colours) && ischar(colours)
            switch colours
                case 'DEC'
                    v_idx = repmat(1:n_pts, [nv 1]);
                    v_idx = v_idx(:);
                    C{i} = abs(rtv(v_idx, :));
                case 'random'
                    slc = rand(1, 3);
                    C{i} = repmat(slc, [n_pts * nv 1]);
            end
        elseif ~isempty(colours) && isnumeric(colours) && size(colours, 1) == nsl
            C{i} = repmat(colours(i, :), [n_pts * nv 1]);
        else                 
            C{i} = ones(n_pts * nv, 1) * i;
        end
    end
        
    F2 = F{1};
    anv = 0;
    for i = 2:nsl
        anv = anv + size(V{i - 1}, 1);
        F2 = [F2; F{i} + anv];
    end
    F = F2; 
    V = cell2mat(V);
    C = cell2mat(C);
    
end

function [v, f] = gen_disc(nf, r)

    % generate a mesh representation of a disc
    alpha = 2 * pi / nf; 
    v = [sin(alpha .* (0:nf - 1)); cos(alpha .* (0:nf - 1)); zeros(1, nf)]';
    v = (v - repmat(mean(v), [size(v, 1) 1])) * r;
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