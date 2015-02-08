% % ps1

%% 1-a
img = imread(fullfile('input', 'ps1-input0.png'));  % already grayscale
%% TODO: Compute edge image img_edges

img_edges = edge(img,'canny');

imwrite(img_edges, fullfile('output', 'ps1-1-a-1.png'));  % save as output/ps1-1-a-1.png

%% 2-a
[H, theta, rho] = hough_lines_acc(img_edges);  % defined in hough_lines_acc.m
%% TODO: Plot/show accumulator array H, save as output/ps1-2-a-1.png
imshow(imadjust(mat2gray(H)),'XData',theta,'YData',rho,...
      'InitialMagnification','fit');
axis on, axis normal, hold on;
colormap(hot);
xlabel('\theta')
ylabel('\rho')
title({'ps1-2-a-1.png','Hough transform of ps1-1-a-1.png'})
hold off;
%pause;

%imwrite(imadjust(mat2gray(H)), fullfile('output', 'ps1-2-a-1.png'));

%% 2-b
peaks = hough_peaks(H, 10);  % defined in hough_peaks.m
%% TODO: Highlight peak locations on accumulator array, save as output/ps1-2-b-1.png
figure;
imshow(imadjust(mat2gray(H)),'XData',theta,'YData',rho,...
      'InitialMagnification','fit');
axis on, axis normal, hold on;
colormap(hot);
xlabel('\theta')
ylabel('\rho')
plot(theta(peaks(:,2)),rho(peaks(:,1)),'bo')
legend('peaks','Location','NorthEast')
title({'ps1-2-b-1.png','Hough transform with peaks'})

%% TODO: Rest of your code here

%% 2-c
hough_lines_draw(img,'ps1-2-c-1.png',peaks,rho,theta);

%% 3
img_noisy = imread(fullfile('input', 'ps1-input0-noise.png'));

%% 3-a smoothing
f = fspecial('gaussian',[3 3],2);

img_noisy_smoothed = imfilter(img_noisy,f);
imwrite(img_noisy_smoothed, fullfile('output', 'ps1-3-a-1.png'));

%% 3-b edge detection
img_edges_noisy = edge(img_noisy,'canny',[0.2 0.55],2);
imwrite(img_edges_noisy, fullfile('output', 'ps1-3-b-1.png'));
img_edges_noisy_smoothed = edge(img_noisy_smoothed,'canny',[0.2 0.55],2);
imwrite(img_edges_noisy_smoothed, fullfile('output', 'ps1-3-b-2.png'));

%% 3-c hough with noise
[H_ns, theta_ns, rho_ns] = hough_lines_acc(img_edges_noisy_smoothed);
peaks_ns = hough_peaks(H_ns, 10, 'Threshold', 0.3*max(H_ns(:)));
figure;
imshow(imadjust(mat2gray(H_ns)),'XData',theta_ns,'YData',rho_ns,...
      'InitialMagnification','fit');
axis on, axis normal, hold on;
colormap(hot);
xlabel('\theta')
ylabel('\rho')
plot(theta_ns(peaks_ns(:,2)),rho_ns(peaks_ns(:,1)),'bo')
legend('peaks','Location','NorthEast')
title({'ps1-3-c-1.png','Hough transform of noisy image with peaks'})

hough_lines_draw(img_noisy,'ps1-3-c-2.png',peaks_ns,rho_ns,theta_ns);

%% 4 coins and pens
img_pens = rgb2gray(imread(fullfile('input', 'ps1-input1.png')));

%% 4-a smoothed
f = fspecial('gaussian',[3 3],2);
img_pens_smoothed = imfilter(img_pens,f);
imwrite(img_pens_smoothed, fullfile('output', 'ps1-4-a-1.png'));

%% 4-b edge detection
img_edges_ps = edge(img_pens_smoothed,'canny',[0.2 0.55],2);
imwrite(img_edges_ps, fullfile('output', 'ps1-4-b-1.png'));

%% 4-c hough for pens
[H_ps, theta_ps, rho_ps] = hough_lines_acc(img_edges_ps);
peaks_ps = hough_peaks(H_ps, 10, 'Threshold', 0.3*max(H_ps(:)),'NHoodSize',floor(size(H) / 110.0) * 2 + 1);
figure;
imshow(imadjust(mat2gray(H_ps)),'XData',theta_ps,'YData',rho_ps,...
      'InitialMagnification','fit');
axis on, axis normal, hold on;
colormap(hot);
xlabel('\theta')
ylabel('\rho')
plot(theta_ps(peaks_ps(:,2)),rho_ps(peaks_ps(:,1)),'bo')
legend('peaks','Location','NorthEast')
title({'ps1-4-c-1.png','Hough transform with peaks'})

hough_lines_draw(img_pens,'ps1-4-c-2.png',peaks_ps,rho_ps,theta_ps);

%% 5 circles
img_pens = rgb2gray(imread(fullfile('input', 'ps1-input1.png')));

%% 5-a
f = fspecial('gaussian',[3 3],2);
img_pens_smoothed = imfilter(img_pens,f);
imwrite(img_pens_smoothed, fullfile('output', 'ps1-5-a-1.png'));

img_edges_ps = edge(img_pens_smoothed,'canny',[0.1 0.55],1.5);
imwrite(img_edges_ps, fullfile('output', 'ps1-5-a-2.png'));

[centers,radii] = find_circles(img_edges_ps,[20 20]);
hough_circles_draw(img_pens,'ps1-5-a-3.png',centers,radii);

%% 5-b
[centers,radii] = find_circles(img_edges_ps,[20 50]);
hough_circles_draw(img_pens,'ps1-5-b-1.png',centers,radii);

%% display circle hough transform for radius=20
radius = 20;
H_circ = hough_circles_acc(img_edges_ps, radius);
peaks_circ = hough_peaks(H_circ, 10);
figure;
imshow(imadjust(mat2gray(H_circ)),...
      'InitialMagnification','fit');
axis on, axis normal, hold on;
colormap(hot);
xlabel('x of center of circle')
ylabel('y of center of circle')
plot(peaks_circ(:,2),peaks_circ(:,1),'bo')
legend('peaks','Location','NorthEast')

% 6
img_input2 = rgb2gray(imread(fullfile('input', 'ps1-input2.png')));
f = fspecial('gaussian',[3 3],2);
img_input2_smoothed = imfilter(img_input2,f);
img_edges_input2 = edge(img_input2_smoothed,'canny',[0.2 0.6],1.5);

% figure; imshow(img_edges_input2)

[H2, theta2, rho2] = hough_lines_acc(img_edges_input2);
peaks2 = hough_peaks(H2,10,'Threshold',0.3*max(H2(:)),'NHoodSize',floor(size(H2) / 200.0) * 2 + 1);

figure;
imshow(imadjust(mat2gray(H2)),'XData',theta2,'YData',rho2,...
      'InitialMagnification','fit');
axis on, axis normal, hold on;
colormap(hot);
xlabel('\theta')
ylabel('\rho')
plot(theta2(peaks2(:,2)),rho2(peaks2(:,1)),'bo')
legend('peaks','Location','NorthEast')
%title({'ps1-4-c-1.png','Hough transform with peaks'})

hough_lines_draw(img_input2,'ps1-6-a-1.png',peaks2,rho2,theta2);

peaks2_pens = [];
pkthetas = theta2(peaks2(:,2));
for n=1:length(pkthetas)
    pktheta = pkthetas(n);
    for m=n+1:length(pkthetas)
        % if parallel and not horizontal or vertical
        if abs(pktheta - pkthetas(m)) < 1.5 && pktheta > 5
            peaks2_pens = [peaks2_pens; peaks2([n m],:)];
        end
    end
end

hough_lines_draw(img_input2,'ps1-6-c-1.png',peaks2_pens,rho2,theta2);

%% 7-a
[centers2,radii2] = find_circles(img_edges_input2,[20 50]);
hough_circles_draw(img_input2,'ps1-7-a-1.png',centers2,radii2);

%% 8-a
img_distort = rgb2gray(imread(fullfile('input', 'ps1-input3.png')));
img_edge_distort = edge(img_distort,'canny',[0.05 0.4],1);

%figure; imshow(img_edge_distort)

[HD,thetaD,rhoD] = hough_lines_acc(img_edge_distort);
peaksDL = hough_peaks(HD,8,'Threshold',0.27*max(H2(:)),'NHoodSize',floor(size(H2) / 200.0) * 2 + 1);
HDlineImg = hough_lines_draw(img_distort,'',peaksDL,rhoD,thetaD);

[centersHD,radiiHD] = find_circles(img_edge_distort,[20 50]);
HDcircImg = hough_circles_draw(img_distort,'',centersHD,radiiHD);

combImg = img_distort;
combImg(HDlineImg(:,:,2)==255 | HDcircImg(:,:,2) == 255) = 255;
HDimg = HDlineImg;
HDimg(:,:,2) = combImg;
imwrite(HDimg, fullfile('output', 'ps1-8-a-1.png'));

%% 8-c
distortion_angle = 50;

%display hough transform for ellipses
radius=50;
H_ellipse = hough_ellipses_acc(img_edge_distort, radius, distortion_angle);
peaks_ellipse = hough_peaks(H_ellipse, 10);
figure;
imshow(imadjust(mat2gray(H_ellipse)),...
      'InitialMagnification','fit');
axis on, axis normal, hold on;
colormap(hot);
xlabel('x of center of ellipse')
ylabel('y of center of ellipse')
plot(peaks_ellipse(:,2),peaks_ellipse(:,1),'bo')
title({'Hough transform: Ellipses','distortion angle \phi=50, major axis=50'})
legend('peaks','Location','SouthEast')

[centersE,radiiE] = find_ellipses(img_edge_distort,[20 50],distortion_angle);
HDellipsImg = hough_ellipses_draw(img_distort,'',centersE,radiiE,distortion_angle);

combImg = img_distort;
combImg(HDlineImg(:,:,2)==255 | HDellipsImg(:,:,2) == 255) = 255;
HDEimg = HDlineImg;
HDEimg(:,:,2) = combImg;
imwrite(HDEimg, fullfile('output', 'ps1-8-c-1.png'));
