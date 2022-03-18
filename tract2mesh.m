function [V, F, DEC] = bundle2mesh
    
    %generate_polygonal_cylinder

    r = 0.1;        % radius
    nv = 6;         % number of vertices at cross-section
    
    tracks = simulate_streamline_bundle([], 1);
    [V, F, DEC] = deal([]);    
    anv = 0;        % current number of vertices to add
    
    % generate cross-sectional disc    
    [dv, df] = gen_disc(nv, r);  
    
    % generate cylindroids
    for i = 1:numel(tracks)
       
        sl = tracks{i};
        n_pts = size(sl, 1);
        
        % gradient vectors - used for rotation and colour-coding
        rtv = norm_vec([gradient(sl(:, 1)) gradient(sl(:, 2)) gradient(sl(:, 3))]);
        norot = all(abs(rtv) == repmat([0 0 1], [n_pts 1]), 2);
        
        % rotation matrices for each cross-section corresponding to each streamline point
        cross_v = cross(repmat([0 0 1], [n_pts, 1]), rtv, 2);
        dot_v = dot(repmat([0 0 1], [n_pts, 1]), rtv, 2);
        
        % position the first disc
        if ~norot(1)            
            R = rot_mtx(cross_v(1, :), dot_v(1, :), sl(1, :));
            v1 = (R * [dv ones(nv, 1)]')';
            
        else
            v1 = dv;
        end
        V = [V; v1(:, 1:3)];
        F = [F; df + anv];
        DEC = [DEC; repmat(rtv(1, :), [nv 1])];
       
        for j = 2:n_pts
            
            % position the vertices of the next disc
            if ~norot(j)
                R = rot_mtx(cross_v(j, :), dot_v(j, :), sl(j, :));
                v1 = (R * [dv ones(nv, 1)]')';
            else
                v1 = dv;
            end 
            
            V = [V; v1(:, 1:3)];
            DEC = [DEC; repmat(rtv(j, :), [nv 1])];
            
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
        
    end
    
    DEC = abs(DEC);
    
end

function [v, f] = gen_disc(nf, r)

    alpha = 2 * pi / nf; 
    v = [sin(alpha .* (0:nf - 1)); cos(alpha .* (0:nf - 1)); zeros(1, nf)]';
    v = (v - mean(v)) * r;
    f = [ones(1, nf - 2); 2:nf - 1; 3:nf]';        

end

function R = rot_mtx(cv, dv, O)

    ssc = [0 -cv(3) cv(2); cv(3) 0 -cv(1); -cv(2) cv(1) 0];
    R = eye(3) + ssc + ssc ^ 2 * (1 - dv) / sqrt(sum(cv .^ 2)) ^ 2;
    R = [R(1, :) O(1); R(2, :) O(2); R(3, :) O(3); 0 0 0 1];    

end

function v = norm_vec(v)

    v = v ./ repmat(sqrt(sum(v .^ 2, 2)), [1 3]);

end