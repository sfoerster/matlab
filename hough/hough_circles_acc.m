function H = hough_circles_acc(BW, radius)
    % Compute Hough accumulator array for finding circles.
    %
    % BW: Binary (black and white) image containing edge pixels
    % radius: Radius of circles to look for, in pixels

    % TODO: Your code here
    
    [M,N] = size(BW);
    H = zeros(M,N);
    
    numpts = ceil(2*pi*radius);
    theta = linspace(0,2*pi,numpts);
    
    y=0;
    for m=M:-1:1
        y=y+1;
        for n=1:N
            x=n;
            
            if BW(m,n)
                
                allx = round(radius.*cos(theta))+n;
                ally = round(radius.*sin(theta))+m;
                
                delidx = allx > N | allx < 1;
                delidx = delidx | ally > M | ally < 1;
                
                allx(delidx) = [];
                ally(delidx) = [];
                
                
                setidx = sub2ind(size(H), ally, allx);
                H(setidx) = H(setidx) + 1;
                
            end
            
        end
    end
end
