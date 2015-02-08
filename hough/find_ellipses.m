function [centers, radii] = find_ellipses(BW, radius_range, distortion_angle)
    % Find circles in given radius range using Hough transform.
    %
    % BW: Binary (black and white) image containing edge pixels
    % radius_range: Range of circle radii [min max] to look for, in pixels

    % TODO: Your code here
    
    centers = [];
    radii = [];
%     for r=radius_range(2):-1:radius_range(1)
    for r=radius_range(1):radius_range(2)
        disp(['[find_ellipses] Radius: ' num2str(r)])
        
        H = hough_ellipses_acc(BW,r,distortion_angle);
        peaks = hough_peaks(H,10,'Threshold',50);
        
        if ~isempty(peaks)
            for m=1:size(peaks,1)
                addPeak = true;
                for k=1:size(centers,1)
                    d = sqrt( (centers(k,1)-peaks(m,1)).^2 + (centers(k,2)-peaks(m,2)).^2 );
                    %disp(['r: ' num2str(r) ' d: ' num2str(d)]) 
                    if d < radii(k)+r
                        addPeak = false;
                    end
                end
                
                if addPeak
                    centers = [centers; peaks(m,:)];
                    radii = [radii; r];
                end
            end
            
%             centers = [centers; peaks];
%             radii = [radii; repmat(r,size(peaks,1),1)];
        end
    end
end

% function overlap = find_circle_overlap(center1,radius1,center2,radius2)
% 
%     d = sqrt( (center1(1)-center2(1)).^2 + (center1(2)-center2(2)).^2 );
%     
%     overlap = (radius1.^2) .* acos( (d.^2 + radius1.^2 - radius2.^2) ./ (2.*d.*radius1)) + ...
%         (radius2.^2) .* acos( (d.^2 + radius2.^2 - radius1.^2) ./ (2.*d.*radius2)) - ...
%         0.5 .* sqrt( (-d+radius1+radius2).*(d+radius1-radius2).* ...
%         (d-radius1+radius2).*(d+radius1+radius2));
% 
% end
