function C = convert_colourmap(C, map, varargin)

    % C:        some scalar values per vertex, Nx1 or Nx3
    % map:      desired colourmap, only applied to Nx1 inputs
    %           - provided as either an Mx3 RGB matrix or as sa name of a pre-defined map 
    %             (type: doc colormap, see options under colormap name)
    %           - if not provided, C is just tidied up for output    
    % varargin: min-max scale
    %
    % examples: C = convert_colourmap(C, []); - will just make a 3-column matrix if needed
    %           C = convert_colourmap(C, [], [0 1]) - as above but will rescale values to [0 1]
    %           C = convert_colourmap(C, [], [min(C(:)) max(C(:))]) - as above but will rescale values to min-max of C
    %           C = convert_colourmap(C, 'jet'); - will allocate jet colourmap values to C scaled at min-max of C
    %           C = convert_colourmap(C, 'jet', [0 1]) - as above but scaled at [0 1]    
    %           M = rand(16, 3); C = convert_colourmap(C, M); - will generate a colourmap of 16 random colours then allocate them to C based on C values
    
    % sort if cell
    isc = false;
    asl = [];
    if iscell(C)
        isc = true;
        sll = cellfun(@length, C);
        asl = [0 cumsum(sll)];
        C = cell2mat(C(:));
    end
    
    % set up scaling
    if nargin > 2 && ~isempty(varargin{1})
        scl = varargin{1};
    else
        scl = [min(C(:)) max(C(:))];
    end
    
    % in case no mapping is needed
    if isempty(map) || size(C, 2) == 3
        
        % rescale if requested
        if nargin > 2 && ~isempty(varargin{1})
            scl = varargin{1};
            C = (C - scl(1)) / scl(2);
        end
        
        % expand to three columns as required
        if size(C, 2) == 1
            C = repmat(C, [1 3]);
        end
        
        check_range(C);    
        C = back2cell(C, asl, isc);
        return
        
    end
        
    % get the pre-defined map if needed
    if ischar(map)
        set(0,'DefaultFigureVisible','off'); % suppresses the figure that pops up
        map = colormap(map);        
        map = map / max(map(:)); % place on a [0 1] scale
        set(0,'DefaultFigureVisible','on');
    end
    
    % rescale
    C = (C - scl(1)) / (scl(2) - scl(1));
    
    % index 
    bins = linspace(0, 1, size(map, 1));
    [~, binned_c] = histc(C, bins);
    C = map(binned_c, :);
    
    check_range(C);
    C = back2cell(C, asl, isc);
    
end

function check_range(C)

    range = C > 1 | C < 0;
    if any(range(:))
        warning('C is out of [0 1] range')
    end

end

function C2 = back2cell(C, asl, isc)

    if ~isc
        C2 = C;
        return
    end
    
    nsl = numel(asl) - 1;
    C2 = cell(nsl, 1);
    for i = 1:nsl
        C2{i} = C(asl(i)+1:asl(i+1), :);
    end
end