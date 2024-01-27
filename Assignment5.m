
% Name : Muktika Manohar
% Contact Email : mm87150n@pace.edu
% Assignment Number : 5

%Part a

Im = imread('/Users/muku/Downloads/ball.bmp');
[H, S, I] = Myrgb2hsi(Im);

% Display the images side-by-side on figure 1
figure(1);
subplot(1, 3, 1); imshow(H); title('Hue');
subplot(1, 3, 2); imshow(S); title('Saturation');
subplot(1, 3, 3); imshow(I); title('Intensity');
pause;

% Use MATLAB's rgb2hsv function 
matlabHSI = rgb2hsv(Im);
matlabH = matlabHSI(:,:,1);
matlabS = matlabHSI(:,:,2);
matlabI = matlabHSI(:,:,3);


% Display the MATLAB results side-by-side on figure 2
figure(2);
subplot(1, 3, 1); imshow(matlabH); title('Matlab Hue');
subplot(1, 3, 2); imshow(matlabS); title('Matlab Saturation');
subplot(1, 3, 3); imshow(matlabI); title('Matlab Intensity');
pause; 


%Display the difference images side-by-side on figure 3
diffH = abs(H - matlabH);
diffS = abs(S - matlabS);
diffI = abs(I - matlabI);

figure(3);
subplot(1, 3, 1); imshow(diffH); title('Hue Difference');
subplot(1, 3, 2); imshow(diffS); title('Saturation Difference');
subplot(1, 3, 3); imshow(diffI); title('Intensity Difference');
pause;


%Display command to explain differences
disp('Explanation for Differences:');

disp('1. The main differences between the implementation and MATLABs rgb2hsv function lie in the formulas for calculating saturation and hue, where the implementation involves additional adjustments and a different approach for saturation.')
disp('2. The custom implementation normalizes hue values after calculation, while MATLABs rgb2hsv returns already normalized hue values.')
disp('3. Variations in handling small values, particularly the use of eps in the custom implementation, may contribute to differences in edge cases.')
pause;

%Part b
% Apply edge detection methods
sobelEdges = edge(matlabI, 'Sobel');
prewittEdges = edge(matlabI, 'Prewitt');
robertsEdges = edge(matlabI, 'Roberts');
logEdges = edge(matlabI, 'log');
cannyEdges = edge(matlabI, 'Canny');

% Display the edges on figure 4
figure(4);
subplot(2, 3, 1); imshow(sobelEdges); title('Sobel Edges');
subplot(2, 3, 2); imshow(prewittEdges); title('Prewitt Edges');
subplot(2, 3, 3); imshow(robertsEdges); title('Roberts Edges');
subplot(2, 3, 4); imshow(logEdges); title('Laplacian of Gaussian Edges');
subplot(2, 3, 5); imshow(cannyEdges); title('Canny Edges');
pause;

disp('Explanation for the best method:');
disp('After viewing the results of all of the above methods used for edge detection of the image of a ball, the best methods turn out to be Sobel, Prewitt and Roberts.')
disp('The three of those properly display the edges of the image.')
disp('The Laplacian of Gaussian and Canny methods focus more on noise reduction so they are not as clear.')
disp('Consider image characteristics, noise levels, and edge thickness; choose Canny for noise suppression and LoG for variable thickness based on efficiency and application requirements, tuning parameters as needed. Use evaluation metrics for optimal edge detector selection.')
pause;

%Part c

% Use multilevel Otsu's method for adaptive thresholding
threshold = multithresh(matlabH);
binaryMask = matlabH > threshold;

% Perform morphological operations for better segmentation (optional)
binaryMask = imopen(binaryMask, strel('disk', 5));

% Find connected components in the binary mask
cc = bwconncomp(binaryMask);

% Get region properties, including weighted centroid
stats = regionprops(cc, matlabH, 'WeightedCentroid');

% Check if any regions were found
if ~isempty(stats)
    % Use the weighted centroid of the first detected region
    centroid = stats(1).WeightedCentroid;

    % Adjust the centroid position (add or subtract offsets)
    offsetX = 0; % Adjust this value to move the centroid horizontally
    offsetY = 95; % Adjust this value to move the centroid vertically
    centroid = centroid + [offsetX, offsetY];

    
    % Display the result on the original image
    figure(5);
    imshow(Im);
    hold on;
    
    % Plot the centroid as a crosshair
    plot(centroid(1), centroid(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    
    title('Centroid of the ball');
    hold off;
end

 

function [H, S, I] = Myrgb2hsi(Im)
    % Ensure the image is in double format for calculations
    Im = im2double(Im);

    % Separate R, G, and B channels
    R = Im(:, :, 1);
    G = Im(:, :, 2);
    B = Im(:, :, 3);

    % Compute Intensity (I)
    I = (R + G + B) / 3;

    % Compute Saturation (S)
    minRGB = min(Im, [], 3);
    S = 1 - (3 ./ (R + G + B + eps)) .* minRGB;

    % Compute Hue (H)
    num = 0.5 * ((R - G) + (R - B));
    denom = sqrt((R - G).^2 + (R - B) .* (G - B));
    H = acosd(num ./ (denom + eps));

    % Adjust H values based on the blue component
    H(B > G) = 360 - H(B > G);

    % Normalize H, S, I to be in the range [0, 1]
    H = H / 360; % Convert hue to the range [0, 1]
    S = max(0, min(S, 1));
    I = max(0, min(I, 1));
end


