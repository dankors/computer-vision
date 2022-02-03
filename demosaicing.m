% Part (a)
% BE SURE TO CHANGE TO PROPER PATH
original = imread('/Users/danielkorsunsky/Desktop/McGill/Fall 2020/COMP 558/Assignment 1/Q1 (Demosaicing)/rainbow.jpg');
imshow(original)

[rows, columns, ~] = size(original)
redChannel = zeros(rows, columns, 'uint8');
blueChannel = zeros(rows, columns, 'uint8');
greenChannel = zeros(rows, columns, 'uint8');

% Sample image simulating Bayer pattern
redChannel = original(:, :, 1)
redChannel(1:2:end, :) = 0
redChannel(:, 2:2:end) = 0

greenChannel = original(:, :, 2)
greenChannel(1:2:end, 2:2:end) = 0
greenChannel(2:2:end, 1:2:end) = 0

blueChannel = original(:, :, 3)
blueChannel(:, 1:2:end) = 0
blueChannel(2:2:end, :) = 0

imshow(redChannel)
imshow(greenChannel)
imshow(blueChannel)


% Part (b)
rbFilter = [.25 .5 .25; .5 1 .5; .25 .5 .25]
gFilter = [0 .25 0; .25 1 .25;  0 .25 0]

redFixed = uint8(conv2(redChannel, rbFilter, 'same'))
blueFixed = uint8(conv2(blueChannel, rbFilter, 'same'))
greenFixed = uint8(conv2(greenChannel, gFilter, 'same'))
imshow(redFixed)
imshow(greenFixed)
imshow(blueFixed)

% Combine the 3 interpolated channels into an RGB image
rgbImage = cat(3, redFixed, greenFixed, blueFixed)
imshow(rgbImage)


% Part (c): same thing but with Ukraine image
