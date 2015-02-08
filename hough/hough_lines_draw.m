function outimg = hough_lines_draw(img, outfile, peaks, rho, theta)
    % Draw lines found in an image using Hough transform.
    %
    % img: Image on top of which to draw lines
    % outfile: Output image filename to save plot as
    % peaks: Qx2 matrix containing row, column indices of the Q peaks found in accumulator
    % rho: Vector of rho values, in pixels
    % theta: Vector of theta values, in degrees

    % TODO: Your code here

    outimg(:,:,1) = img;
    outimg(:,:,2) = img;
    outimg(:,:,3) = img;
    
    lineImg = flipud(img);

    for k=1:size(peaks,1)
        
        thisrho = rho(peaks(k,1));
        thistheta = theta(peaks(k,2));

        [M, N] = size(img);
        %lineImg = nan(M,N);

        xyresthresh = 1.5.*sqrt(1/pi);

        m = M:-1:1;
        n = 1:N;
        y = (1:M)';
        x = 1:N;

        [x,y] = meshgrid(x,y);

        rhoTrue = x.*cosd(thistheta) + y.*sind(thistheta);
        %[yout,xout] = find(abs(rhoTrue - thisrho) < xyresthresh);
        lineImg(abs(rhoTrue - thisrho) < xyresthresh) = 255;

        %colorout = uint8(zeros(M,N,3));
        %colorout(:,:,2) = uint8(lineImg);

    end
    
    outimg(:,:,2) = flipud(lineImg);
    
    if ~isempty(outfile)
        imwrite(outimg, fullfile('output', outfile));
    end

end
