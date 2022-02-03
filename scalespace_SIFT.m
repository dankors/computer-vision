% Remember to check filepath for images! (shapes.jpg and octopus.jpg)

image = rgb2gray(imread('/Users/danielkorsunsky/Desktop/McGill/Fall 2020/COMP 558/Assignment 2/octopus.jpg'));
%image = rgb2gray(imread('/Users/danielkorsunsky/Desktop/McGill/Fall 2020/COMP 558/Assignment 2/shapes.jpg'));
tau = 4; %4 for octopus, 10 for shapes

% For Q7: Scale invariance (turn on and off with comment)
%image = imresize(image, .5);


% Q1: Gaussian scale space

[rows, columns] = size(image);
gss = zeros(rows, columns, 17);

for s = 0:16
    sigma = 2 ^ (s/4);
    gss(:, :, s+1) = imgaussfilt(image, sigma);
end

for slice = 2:17
    subplot(4,4,slice-1)
    imagesc(gss(:,:,slice)); colormap("gray");
end


% Q2: Harris-Stevens

hs = cell(1,17);

for level = 1:17
    slice = gss(:,:,level);
    dIdx = conv2(slice,  [1 0 -1]);
    dIdy = conv2(slice,  [1 0 -1]');
    dIdx = dIdx(:,2:size(dIdx,2)-1); % made derivatives same size
    dIdy = dIdy(2:size(dIdy,1)-1,:);

    Mxx = dIdx.^2;
    Mxy = dIdx .* dIdy;
    Myy = dIdy.^2;
    
    %sigma of ngd is twice the size of sigma of scale space slice
    sigma = 2 ^ ((level-1)/4);
    width = round(10 * sigma);
    weight = fspecial('gaussian', width, 2*sigma);

    % weight of second moment matrix
    Mxx = conv2(Mxx, weight, 'valid');
    Mxy = conv2(Mxy, weight, 'valid');
    Myy = conv2(Myy, weight, 'valid');

    det = Mxx .* Myy - Mxy.^2;
    tr = Mxx + Myy;

    harrisStevens = det - .1 * (tr.^2);
    hs{level} = harrisStevens;
end

close all
for level = 2:17
    subplot(4,4,level-1)
    imagesc(hs{level}); axis off; colormap("jet"); colorbar
end


% Q3: Difference of Gaussians

DOG = zeros(size(gss));

for j = 1:16
   DOG(:,:,j) = gss(:,:,j) - gss(:,:,j+1);
end

close all
for layer = 1:16
    subplot(4,4,layer)
    imagesc(DOG(:,:,layer)); colormap("jet"); colorbar
end


% Q4: SIFT Keypoints

[ySize, xSize] = size(DOG(:,:,1));
keypoints = zeros(1,4);

for layer = 2:15
    sigma = 2 ^ ((layer-1)/4);
    bound = ceil(2*sigma); % at least 2*sigma from image boundary
    
    for x = bound:xSize - bound
        for y = bound:ySize - bound
            % PIQ = point in question
            
            mag = abs(DOG(y, x, layer)); % magnitude of PIQ
            neighborhood = DOG(y-1:y+1, x-1:x+1, layer-1:layer+1); %3x3x3 ngd
            maxes = find(neighborhood == max(max(max(neighborhood))));
            mins = find(neighborhood == min(min(min(neighborhood))));

            % Max/min has to be unique, correspond to PIQ (index = 14 in 
            % 3x3x3 matrix), and have magnitude larger than threshold tau
            if size(maxes, 1) == 1 && maxes(1) == 14 && mag > tau
                keypoints = [keypoints; [x y sigma 0]];
                break
            end
            if size(mins, 1) == 1 && mins(1) == 14 && mag > tau
                keypoints = [keypoints; [x y sigma 0]];
                break
            end
            
        end
    end
end
keypoints(1, :) = []; %remove row of zeros from beginning (needed for init)

% Plot keypoints onto image
close all
imshow(image)
hold on
for k = 1:size(keypoints, 1)
    kp = keypoints(k, :); x = kp(1); y = kp(2); s = kp(3);
    viscircles([x y], s, 'Color', 'red')
end


% Q5: Hessian constraint

r = 10; % as suggested by Lowe (2004)
threshold = ((r + 1)^2) / r;

for level = 2:15
    slice = DOG(:,:,level);
    dIdx = conv2(slice,  [1 0 -1], 'same'); 
    dIdy = conv2(slice,  [1 0 -1]', 'same');
    
    % approximation of Hessian 2nd partial derivatives w/ local diff filter
    Dxx = conv2(dIdx,  [1 0 -1], 'same'); 
    Dxy = conv2(dIdx,  [1 0 -1]', 'same'); 
    Dyy = conv2(dIdy,  [1 0 -1]', 'same');
    
    det = Dxx .* Dyy - Dxy.^2;
    tr = Dxx + Dyy;

    sigma = 2 ^ ((level-1)/4);
    % iterate through keypoints backwards because length of keypoints list
    % decreases with removal 
    for k = size(keypoints, 1):-1:1
        kp = keypoints(k, :);   
        
        if kp(3) == sigma % only check keypoints on this sigma level
            x = kp(1); y = kp(2);
            % remove keypoint if det <= o
            if det(y, x) <= 0
                keypoints(k, :) = [];
                continue
            end
            
            % remove keypoint if ratio > threshold
            ratio = (tr(y,x)^2) / det(y,x);           
            if ratio > threshold
                keypoints(k, :) = [];
            end
        end
    end
end

%close all; imshow(image)
hold on
for k = 1:size(keypoints, 1)
    kp = keypoints(k, :); x = kp(1); y = kp(2); s = kp(3);
    viscircles([x y], s, 'Color', 'cyan')
end


% Q6: SIFT feature dominant orientation

for layer = 2:15
    slice = gss(:,:,layer);
    dIdx = conv2(slice,  [1 0 -1]); dIdy = conv2(slice,  [1 0 -1]');
    dIdx = dIdx(:,2:size(dIdx,2)-1); dIdy = dIdy(2:size(dIdy,1)-1,:);
    
    Imag = sqrt(dIdx .* dIdx + dIdy .* dIdy);
    sigma = 2 ^ ((layer-1)/4);
    Imag = imgaussfilt(Imag, 1.5*sigma);
    Idir = atan2(dIdy, dIdx); Idir = Idir * (180/pi);
    
    for k = 1:size(keypoints, 1)
        orientationHist = zeros(1,36);
        kp = keypoints(k, :);
        
        if kp(3) == sigma
            kp_x = kp(1); kp_y = kp(2); r = floor(2*sigma);
            
            % size of window is 2*sigma
            for y = kp_y-r:kp_y+r
                for x = kp_x-r:kp_x+r
                    theta = Idir(y,x); 
                    if theta < -5 % correct negative degree values
                       theta = theta + 360; 
                    end
                    bin = floor((theta + 15)/10); % convert theta to proper bin
                    orientationHist(bin) = orientationHist(bin) + Imag(y,x);
                end
            end
            
            max_dir = find(orientationHist == max(orientationHist)); % maximum orientation
            keypoints(k, 4) = (max_dir(1) - 1) * 10; % add orientation to keypoint
            other_dir = find(orientationHist >= max(orientationHist) * 0.8);% other orientations >= 80% of max
            if ~isempty(other_dir)
                % create new keypoints with other orientations
                for m = 1:length(other_dir)
                    if m ~= max_dir(1) % (max is already added to keypoints)
                        theta = (other_dir(m) - 1) * 10;
                        keypoints = [keypoints; [kp_x kp_y sigma theta]];
                    end
                end
            end
            
            % Orientation histogram for 2 chosen points
            
            % x-axis points
            orientations = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350];
            
            % chosen point for shapes.jpg
            if (kp_x == 196 && kp_y == 92) && tau == 10 
                close all
                plot(orientations, orientationHist, '*-')
                xlabel('Orientation');
                ylabel('Sum of weighted gradient magnitudes');
                title('Orientation Histogram');
            end
            
            % chosen point for octopus.jpg
            if (kp_x == 223 && kp_y == 87) && tau == 4
                close all
                plot(orientations, orientationHist, '*-')
                xlabel('Orientation');
                ylabel('Sum of weighted gradient magnitudes');
                title('Orientation Histogram');
            end
            
        end
    end
end


% Annotate points with orientations
close all
imagesc(image); colormap(gray);
hold on
for k = 1:size(keypoints, 1)
    kp = keypoints(k, :); x = kp(1); y = kp(2); s = kp(3); t = kp(4);
    viscircles([x y], s, 'Color', 'cyan')
end

for k = 1:size(keypoints, 1)
    kp = keypoints(k, :); x = kp(1); y = kp(2); s = kp(3); t = kp(4);
    x2 = round(x + (2*s * cosd(t)));
    y2 = round(y + (2*s * sind(t)));
    line([x,x2], [y,y2], 'Color', 'r', 'LineWidth', 2);
end
