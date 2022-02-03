I1 = imread('/Users/danielkorsunsky/Desktop/McGill/Fall 2020/COMP 558/Assignment 3/IMAGES/stereo1.jpg');
I2 = imread('/Users/danielkorsunsky/Desktop/McGill/Fall 2020/COMP 558/Assignment 3/IMAGES/stereo2.jpg');
I1bw = rgb2gray(I1); 
I2bw = rgb2gray(I2);

ptsImage1 = detectSURFFeatures(I1bw);
ptsImage2 = detectSURFFeatures(I2bw);

[featuresImage1, validPtsImage1] = extractFeatures(I1bw, ptsImage1);
[featuresImage2, validPtsImage2] = extractFeatures(I2bw, ptsImage2);

indexPairs = matchFeatures(featuresImage1, featuresImage2);
matchedPoints1 = validPtsImage1(indexPairs(:,1));
matchedPoints2 = validPtsImage2(indexPairs(:,2));

[F, inliers] = estimateFundamentalMatrix(matchedPoints1,matchedPoints2,'Method','RANSAC','NumTrials',5000,'DistanceThreshold',0.00001);

numberInliers = length(find(inliers == 1));
inliersI1 = matchedPoints1(inliers,:).Location;
inliersI2 = matchedPoints2(inliers,:).Location;

colors = ["red","green","#2460ED","cyan","magenta","yellow","#17962E","#FFA500"];

% The following matrix will be used to find the intersection of two of the
% epipolar lines, which is equivalent to the epipole of image 1.
e1_lines = zeros(2,3);

figure(1); imshow(I1bw); hold on;
for i = 1:numberInliers
    x = inliersI1(i,1);
    y = inliersI1(i,2);
    c = mod(i-1, 8) + 1; % index of color
    plot(x, y, 'rs', 'MarkerSize', 20, 'LineWidth', 3, 'Color', colors(c));
    p1 = [x y 1]';
    l2 = F * p1;
    m = l2(1,1)/l2(2,1); % slope of line
    x1 = 1;
    x2 = size(I1, 2);
    y1 = y - ((x-x1) * m);
    y2 = y + ((x2-x) * m);
    line([x1 x2], [y1 y2], 'LineWidth', 2, 'LineStyle', '--', 'Color', colors(c));
    
    % Find standard form of first two epipolar lines for system of 2 equations
    if (i == 1) || (i == 2)
        coeff = polyfit([x1 x2], [y1 y2], 1);
        e1_lines(i,:) = [(-1 * coeff(1)) 1 coeff(2)];
    end
end
e1 = rref(e1_lines); % solution is coordinates of epipole
str = sprintf('epipole in image 1 is (x,y) = (%d,%d)', round(e1(1,3)), round(e1(2,3)));
title(str);
saveas(gcf, "epipolarlines1.png", 'png');
close all;

e2_lines = zeros(2,3);

figure(2); imshow(I2bw); hold on; title('image2');
for i = 1:numberInliers
    x = inliersI2(i,1);
    y = inliersI2(i,2);
    c = mod(i-1, 8) + 1;
    plot(x, y, 'rs', 'MarkerSize', 20, 'LineWidth', 2, 'Color', colors(c));
    p2 = [x y 1];
    l1 = p2 * F;
    m = l1(1,1)/l1(1,2);
    x1 = 1;
    x2 = size(I2, 2);
    y1 = y - ((x-x1) * m);
    y2 = y + ((x2-x) * m);
    line([x1 x2], [y1 y2], 'LineWidth', 2, 'LineStyle', '--', 'Color', colors(c));
    
    if (i == 1) || (i == 2)
        coeff = polyfit([x1 x2], [y1 y2], 1);
        e2_lines(i,:) = [(-1 * coeff(1)) 1 coeff(2)];
    end
end
e2 = rref(e2_lines);
str = sprintf('epipole in image 2 is (x,y) = (%d,%d)', round(e2(1,3)), round(e2(2,3)));
title(str);
saveas(gcf, "epipolarlines2.png", 'png');


% Disparity Map

[T1, T2] = estimateUncalibratedRectification(F,matchedPoints1,matchedPoints2,size(I1));
[I1Rect,I2Rect] = rectifyStereoImages(I1,I2,T1,T2);
I1RectBW = rgb2gray(I1Rect);
I2RectBW = rgb2gray(I2Rect);

disparityRange = [-48 80];
disparityMap = disparitySGM(I1RectBW,I2RectBW,'DisparityRange',disparityRange,'UniquenessThreshold',20);

hold off; figure(3); 
imshow(disparityMap,disparityRange);
title('Disparity Map'); colormap parula; colorbar
saveas(gcf, "disparitymap.png", 'png');