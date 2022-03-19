function tracks = simulate_streamline_bundle(nsl, step_size)

    % number of streamlines
    if isempty(nsl)
        nsl = randi(10 ^ 3);
    end
    
    % step size    
    if isempty(step_size)
        step_size = 0.5;
    end    
    
    rad = 3;        % central radius
    md = 10;        % maximum dimension
    n_pts = 1000;   % number of points in each curve
    
    %% generate centroid
    
    x = linspace(-md/2, md/2, n_pts);
    y = gen_rand_poly(md, n_pts);
    z = gen_rand_poly(md, n_pts);
    c = [x' y' z'];
    
    % middle point
    cm = median(c);
    [~, cc] = min(sum((c - repmat(cm, [n_pts 1])) .^ 2, 2));   
    
    % tangent
    ct = c(cc + 1, :) - c(cc, :);
    
    
    %% generate streamlines
    
    % normal vector in random direction perpendicular to tangent
    dirvec = rand(nsl, 3);    
    ctr = repmat(ct, [nsl 1]);
    dirvec = dirvec - ctr .* dot(dirvec, ctr);
    dirvec = dirvec ./ sqrt(sum(dirvec .^ 2, 2));  
    
    % pre-generate degree of distortion along each point
    dist = linspace(-1, 1, n_pts)';
    
    tracks = cell(nsl, 1);
    
    for i = 1:nsl
        
        % copy centroid and shift perpendicular to its middle by up to rad
        tracks{i} = c + repmat(dirvec(i, :) * rand * rad, [n_pts 1]);
        
        % add distortion
        tracks{i} = tracks{i} + repmat(dist, [1 3]) .* 10 .* rand;
        
        % resample to match the desired step size
        sl_len = sqrt(sum(diff(tracks{i}) .^ 2, 2));
        sl_bins = [0; cumsum(sl_len)];        
        
        rn_pts = floor(sum(sl_len) / step_size);
        rcum_len = (0:step_size:(rn_pts - 1) * step_size)';
        
        [~, pt_loc] = histc(rcum_len, sl_bins);        
        
        far_down = (rcum_len - sl_bins(pt_loc)) ./ sl_len(pt_loc);
        tracks{i} = tracks{i}(pt_loc, :) + (tracks{i}(pt_loc + 1, :) - tracks{i}(pt_loc, :)) .* repmat(far_down, [1 3]);
        
    end    

end

function p = gen_rand_poly(sz, n_pts)

    % generates a random polynomial to create a curved shape
    poly_ord = randi([2 10]);
    poly_coef = rand(poly_ord, 1) .* (randi(10, poly_ord, 1) - 1) .* sign(rand(poly_ord, 1) - 0.5);
    
    p = polyval(poly_coef, linspace(-sz/2, sz/2, n_pts));
    p = p - min(p);
    p = p / max(p) * sz;
    
end