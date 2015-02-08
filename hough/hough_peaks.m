function peaks = hough_peaks(H, varargin)
    % Find peaks in a Hough accumulator array.
    %
    % Threshold (optional): Threshold at which values of H are considered to be peaks
    % NHoodSize (optional): Size of the suppression neighborhood, [M N]
    %
    % Please see the Matlab documentation for houghpeaks():
    % http://www.mathworks.com/help/images/ref/houghpeaks.html
    % Your code should imitate the matlab implementation.

    %% Parse input arguments
    p = inputParser;
    addOptional(p, 'numpeaks', 1, @isnumeric);
    addParameter(p, 'Threshold', 0.5 * max(H(:)));
    addParameter(p, 'NHoodSize', floor(size(H) / 100.0) * 2 + 1);  % odd values >= size(H)/50
    parse(p, varargin{:});

    numpeaks = p.Results.numpeaks;
    threshold = p.Results.Threshold;
    nHoodSize = p.Results.NHoodSize;

    % TODO: Your code here

    [M, N] = size(H);
    
    peaks = [];
    for m=1:numpeaks
        Hs = sort(H(:),'descend');
        
        pkval = Hs(1);
        %disp(['peak value: ' num2str(pkval)])
        if pkval >= threshold
            [h,k] = find(H==pkval,1);

            % Suppression
            lowX = max([floor(h-nHoodSize(1)) 1]);
            highX = min([ceil(h+nHoodSize(1)) M]);
            lowY = max([floor(k-nHoodSize(2)) 1]);
            highY = min([ceil(k+nHoodSize(2)) N]);
            H(lowX:highX,lowY:highY) = 0;

            peaks = [peaks; h, k];
        end

    end
    

end
