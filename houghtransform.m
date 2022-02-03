% GIVEN SCRIPT   

im = imread("/Users/danielkorsunsky/Desktop/McGill/Fall 2020/COMP 558/Assignment 1/Q3 (Hough)/hallwaytoo.jpg");
im = imresize(rgb2gray(im), 0.5);

Iedges = edge(im,'canny');  
[~,grad_dir]=imgradient(im);
grad_dir = - grad_dir;
imshow(Iedges)

[row, col] = find(Iedges);
edges = [col, row, zeros(length(row),1), zeros(length(row),1) ];
for k = 1:length(row)
     edges(k,3) = cos(grad_dir(row(k),col(k))/180.0*pi);
     edges(k,4) = sin(grad_dir(row(k),col(k))/180.0*pi);
end


% START OF MY CODE
% Create matrix to store (x, y) votes 
pointVotes = zeros(size(Iedges))

% Iterate over edges computed by script
for e = 1:size(edges, 1)
    
    row = edges(e, :)
    x = row(1); y = row(2); cosine = row(3); sine = row(4)
    slope = cosine/(sine)
    
    % If slope is infinity (sine = 0), vote for all (x, y) on vertical line
    if slope == inf || slope == -inf
        pointVotes(:, x) = pointVotes(y, x) + 1
        continue
    end

    % If slope is more vertical than horizontal
    if slope >= 1 || slope <= -1
        
        % Iterate across all y positions from y to 1
        for yPrime = y:-1:1
            xPrime = x + round(slope * (y - yPrime))
            % Break if x-value outside image boundaries
            if slope < 0 && xPrime < 1
                break
            elseif slope > 0 && xPrime > size(pointVotes, 2)
                break
            end
            % Cast vote for (x, y)
            pointVotes(yPrime, xPrime) = pointVotes(yPrime, xPrime) + 1
        end
        
        % Iterate across all y positions from y+1 to max y value
        for yPrime = y+1:size(pointVotes, 1)
            xPrime = x - round(slope * (yPrime - y))
            if slope < 0 && xPrime > size(pointVotes, 2)
                break
            elseif slope > 0 && xPrime < 1
                break
            end
            pointVotes(yPrime, xPrime) = pointVotes(yPrime, xPrime) + 1
        end
        
    % If slope is more horizontal than vertical
    elseif slope < 1 && slope > -1
        for xPrime = x:-1:1
            yPrime = y + round(slope * (x - xPrime))
            if slope < 0 && yPrime < 1
                break
            elseif slope > 0 && yPrime > size(pointVotes, 1)
                break
            end
            pointVotes(yPrime, xPrime) = pointVotes(yPrime, xPrime) + 1
        end
        
        for xPrime = x+1:size(pointVotes, 2)
            yPrime = y - round(slope * (xPrime - x))
            if slope < 0 && yPrime > size(pointVotes, 1)
                break
            elseif slope > 0 && yPrime < 1
                break
            end
            pointVotes(yPrime, xPrime) = pointVotes(yPrime, xPrime) + 1
        end
        
    end
end


imshow(pointVotes)

% Find (x, y) with most votes cast (average coordinates if more than 1 max)
[yMaxes, xMaxes] = find(pointVotes == max(max(pointVotes)))
yMax = round(mean(yMaxes))
xMax = round(mean(xMaxes))

% Plot vanishing point onto image
imagesc(im)
colormap("gray")
hold on
plot(xMax, yMax, 'ro', 'MarkerSize', 15, 'LineWidth', 2);
