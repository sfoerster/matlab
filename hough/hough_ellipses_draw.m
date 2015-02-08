function outimg = hough_ellipses_draw(img, outfile, centers, radii, distortion_angle)

    outimg(:,:,1) = img;
    outimg(:,:,2) = img;
    outimg(:,:,3) = img;
    
    lineImg = img;
    
    for k=1:size(radii,1)
        
        radius = radii(k);
        theta = linspace(0,2*pi,ceil(2*pi*radius));
        x = round(radius.*cos(theta)+centers(k,2));
        y = round(cosd(distortion_angle).*radius.*sin(theta)+centers(k,1));
        
        delidx = x < 1 | x > size(img,2);
        delidx = delidx | y < 1 | y > size(img,1);
        
        x(delidx) = [];
        y(delidx) = [];
        
        setidx = sub2ind(size(lineImg), y, x);
        lineImg(setidx) = 255;
        
    end
    
    outimg(:,:,2) = lineImg;
    
    if ~isempty(outfile)
        imwrite(outimg, fullfile('output', outfile));
    end

end