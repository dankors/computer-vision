% Note: The script may take a while to run, as the conditions for the
% consensus set are rather high (90% of all matched points, with less than
% a 1e-4 error. If it is taking too long, feel free to loosen the
% constraints (Lines 31 and 35).

image1 = imread('/Users/danielkorsunsky/Desktop/McGill/Fall 2020/COMP 558/Assignment 3/IMAGES/cabinets1.jpg');
image2 = imread('/Users/danielkorsunsky/Desktop/McGill/Fall 2020/COMP 558/Assignment 3/IMAGES/cabinets2.jpg');
image1_bw = rgb2gray(image1); image2_bw = rgb2gray(image2);

ptsImage1 = detectSURFFeatures(image1_bw);
ptsImage2 = detectSURFFeatures(image2_bw);

[featuresImage1, validPtsImage1] = extractFeatures(image1_bw, ptsImage1);
[featuresImage2, validPtsImage2] = extractFeatures(image2_bw, ptsImage2);

indexPairs = matchFeatures(featuresImage1, featuresImage2);

matchedImage1 = validPtsImage1(indexPairs(:,1));
matchedImage2 = validPtsImage2(indexPairs(:,2));

numberPairs = length(matchedImage1);
matchedPoints = zeros(numberPairs, 4);
for i = 1:numberPairs
    x1 = matchedImage1.Location(i);
    y1 = matchedImage1.Location(i + numberPairs);
    x2 = matchedImage2.Location(i);
    y2 = matchedImage2.Location(i + numberPairs);
    matchedPoints(i, :) = [x1 y1 x2 y2];
end
C = zeros(1, 4); % initialize consensus set
error = 0.00001;
H = zeros(9,1);
rounds = 0;

while size(C, 1) < 0.8 * numberPairs % at least 90% of points must be in C
    rounds = rounds + 1;
    C = zeros(1, 4);
    sampleIndices = randi(numberPairs, 1, 4);
    samplePoints = zeros(1, 9);

    for sample = 1:4
        index = sampleIndices(sample);
        row = matrixize(matchedPoints, index); % function is at bottom of script
        samplePoints = [samplePoints; row];
    end
    samplePoints(1, :) = [];

    H = null(samplePoints);
    
    for i = 1:numberPairs
        % exclude sample points from homography fit test
        if (i ~= sampleIndices(1)) && (i ~= sampleIndices(2)) && (i ~= sampleIndices(3)) && (i ~= sampleIndices(4))
            row = matrixize(matchedPoints, i);
            prod = row * H;
            aboveError = find(prod >= error);
            if isempty(aboveError)
                C = [C; matchedPoints(i, :)];
            end
        end
    end
    C(1, :) = [];
end

% Data normalization
meanX1 = sum(C(:, 1)) / length(C);
meanY1 = sum(C(:, 2)) / length(C);
meanX2 = sum(C(:, 3)) / length(C);
meanY2 = sum(C(:, 4)) / length(C);

C_temp = zeros(size(C, 1), 2);
C_temp(:, 1) = (C(:, 1) - meanX1) .^ 2 + (C(:, 2) - meanY1) .^ 2;
C_temp(:, 2) = (C(:, 3) - meanX2) .^ 2 + (C(:, 4) - meanY2) .^ 2;
sigma1 = sqrt(sum(C_temp(:,1)) / (2*length(C_temp)));
sigma2 = sqrt(sum(C_temp(:,2)) / (2*length(C_temp)));

C(:, 1) = (C(:, 1) - meanX1) / sigma1;
C(:, 2) = (C(:, 2) - meanY1) / sigma1;
C(:, 3) = (C(:, 3) - meanX2) / sigma2;
C(:, 4) = (C(:, 4) - meanY2) / sigma2;

M1 = [1/sigma1 0 0; 0 1/sigma1 0; 0 0 1] * [1 0 (-1 * meanX1); 0 1 (-1 * meanY1); 0 0 1];
M2 = [1/sigma2 0 0; 0 1/sigma2 0; 0 0 1] * [1 0 (-1 * meanX2); 0 1 (-1 * meanY2); 0 0 1];

C_matrix = zeros(1,9);
for i = 1:length(C)
    row = matrixize(C, i);
    C_matrix = [C_matrix; row];
end
C_matrix(1, :) = [];

[U, S, V] = svd(C_matrix);
H_norm = reshape(V(:, 9), 3, 3)';
H = inv(M2) * H_norm * M1;

% COMPOSITE IMAGE

% set red image into center of composite image
dim = size(image1);
dim(1, 1:2) = 2 * dim(1, 1:2);
composite = zeros(dim);
x1_start = size(composite, 2) / 4; 
y1_start = size(composite, 1) / 4;
x1_end = x1_start + size(image1, 2) - 1; 
y1_end = y1_start + size(image1, 1) - 1;
composite(y1_start:y1_end, x1_start:x1_end, 1) = image1(:,:,1);

% stitch images together
for x = 1:size(composite, 2)
    for y = 1:size(composite, 1)
        I_p = H * [x y 1]';
        x_p = round(I_p(1) / I_p(3));
        y_p = round(I_p(2) / I_p(3));
        if (x_p >= 1) && (x_p <= size(image2, 2)) && (y_p >= 1) && (y_p <= size(image2, 1))
            composite(y + size(image1, 1)/2, x + size(image1, 2)/2, 2:3) = image2(y_p, x_p, 2:3);
        end
    end 
end

imshow(composite/255);


function row = matrixize(pointsList, index)
    x1 = pointsList(index, 1);
    y1 = pointsList(index, 2);
    x2 = pointsList(index, 3);
    y2 = pointsList(index, 4);
    row = [x1 y1 1 0 0 0 (-x2*x1) (-x2*y1) -x2
           0 0 0 x1 y1 1 (-y2*x1) (-y2*y1) -y2];
end
