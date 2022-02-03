% Part(a): Create image of 20 rectangles

image(1:500, 1:500) = 1
for i = 1:20
    xStart = round(100 + 200 * rand)
    yStart = round(100 + 200 * rand)
    xEnd = xStart + round(rand * 150)
    yEnd = yStart + round(rand * 150)
    image(xStart:xEnd, yStart:yEnd) = rand
end
newImage = image(125:374, 125:374)

imagesc(newImage); colormap("gray")

% Part(b)
filt = fspecial('log', 20, 3)
imagesc(filt)
colorbar

% Part(c)
filteredImage = conv2(newImage, filt, 'same')

% Find indices where intensity differences with neighbors are negative
xIndices = find(filteredImage(:, :) .* circshift(filteredImage(:, :), [0 1]) < 0);   
yIndices = find(filteredImage(:, :) .* circshift(filteredImage(:, :), [1 0]) < 0);   
dIndices = find(filteredImage(:, :) .* circshift(filteredImage(:, :), [1 1]) < 0);   
zeroCrossings = zeros(size(filteredImage))
zeroCrossings(xIndices) = 1
zeroCrossings(yIndices) = 1
zeroCrossings(dIndices) = 1

imagesc(filteredImage)
imagesc(zeroCrossings)


% Part(d)
noise = randn(size(newImage)) * (.1 * std2(newImage))
noisyRectangles = newImage + noise
imagesc(noisyRectangles)

noisyFilteredImage = conv2(noisyRectangles, filt, 'same')

xIndices = find(noisyFilteredImage(:, :) .* circshift(noisyFilteredImage(:, :), [0 1]) < 0);   
yIndices = find(noisyFilteredImage(:, :) .* circshift(noisyFilteredImage(:, :), [1 0]) < 0);   
dIndices = find(noisyFilteredImage(:, :) .* circshift(noisyFilteredImage(:, :), [1 1]) < 0);   
zeroCrossings = zeros(size(noisyFilteredImage))
zeroCrossings(xIndices) = 1
zeroCrossings(yIndices) = 1
zeroCrossings(dIndices) = 1

imagesc(noisyFilteredImage)
imagesc(zeroCrossings)


% Part(e)
filt2 = fspecial('log', 40, 6)
imagesc(filt2)
colorbar

fuzzyFilteredImage = conv2(noisyRectangles, filt2, 'same')

xIndices = find(fuzzyFilteredImage(:, :) .* circshift(fuzzyFilteredImage(:, :), [0 1]) < 0);   
yIndices = find(fuzzyFilteredImage(:, :) .* circshift(fuzzyFilteredImage(:, :), [1 0]) < 0);   
dIndices = find(fuzzyFilteredImage(:, :) .* circshift(fuzzyFilteredImage(:, :), [1 1]) < 0);   
zeroCrossings = zeros(size(fuzzyFilteredImage))
zeroCrossings(xIndices) = 1
zeroCrossings(yIndices) = 1
zeroCrossings(dIndices) = 1

imagesc(fuzzyFilteredImage)
imagesc(zeroCrossings)