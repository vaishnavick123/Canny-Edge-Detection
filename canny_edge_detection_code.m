
image = imread('image_mri.jpg');
sigma=3;
threshold=0.01;
highThresholdRatio = 0.275; 
lowThresholdRatio = 0.25;

im = rgb2gray(image);
im = double(im)/255;
imshow(im); title('Input Image');


X = floor(25/2); 
Y = floor(25/2); 
[Xmtx, Ymtx] = meshgrid(-X:X, -Y:Y);  

%Applying Gaussian Filter
gauss = (1/(2*pi*sigma*sigma))*(exp(-(Xmtx.^2 + Ymtx.^2) / (2*sigma*sigma)));
im = imfilter(im, gauss, 'same');    

% Finite difference Img_dx and Img_dy
H = [-1,0,1];
G_dx = conv2(gauss, H, 'valid'); 
figure, subplot(2,2,1), surf(G_dx); title('Gaussian Derivative X-directon');
img_dx = imfilter(im, G_dx, 'same');
H2 = [-1;0;1];
G_dy = conv2(gauss, H2, 'valid');  
subplot(2,2,2), surf(G_dy), title('Gaussian Derivative Y-directon');
img_dy = imfilter(im, G_dy, 'same');



subplot(2,2,3), imshow(img_dx, []), title('Finite Difference:- x Axis');
subplot(2,2,4), imshow(img_dy, []),  title('Finite Difference:- y Axis');

im_gradient = sqrt(img_dx.^2 + img_dy.^2);
figure, imshow(im_gradient, []),  title('Gradient MAGNITUDE');

 
im_threshold = im_gradient > threshold; 
%figure, imshow(im_threshold, []),  title('Gradient Thresholding');

[row col] = size(im_threshold);
non_max_suppression = zeros(row, col);


% Non max suppression 
for i=2:row-1
    for j=2:col-1
        left_px = im_gradient(i-1,j);
        right_px = im_gradient(i+1,j);
        top_px = im_gradient(i,j-1);
        bottom_px = im_gradient(i,j+1);
        center_px = im_gradient(i,j);
        
        if ((center_px > right_px) && (center_px > left_px))
            non_max_suppression(i,j) = center_px;
        end
        
        if ((center_px > top_px) && (center_px > bottom_px))
            non_max_suppression(i,j) = center_px;
        end
    end
end

figure, imshow(non_max_suppression, []),  title('Non-Max Suppression');

%Double Thresholding
highThreshold = max(max(non_max_suppression))*highThresholdRatio;
lowThreshold = highThreshold*lowThresholdRatio;
strongEdgesRow = zeros(1,row*col); 
strongEdgesCol = zeros(1,row*col); 
weakEdgesRow = zeros(1,row*col);  
weakEdgesCol = zeros(1,row*col);  
strongIndex = 1;
weakIndex = 1;
for i=2:row-1
    for j=2:col-1
        if non_max_suppression(i,j) > highThreshold
            non_max_suppression(i,j) = 1;
            strongEdgesRow(strongIndex) = i;
            strongEdgesCol(strongIndex) = j;
            strongIndex = strongIndex + 1;
        elseif non_max_suppression(i,j) < lowThreshold
            non_max_suppression(i,j) = 0;
        else
            weakEdgesRow(weakIndex) = i;
            weakEdgesCol(weakIndex) = j;
            weakIndex = weakIndex + 1;
        end
    end
end
figure; imshow(non_max_suppression);
title('Double Threshold'); 
       
% Edge tracking by hysteresis
set(0,'RecursionLimit',10000)
for i=1:strongIndex-1
    non_max_suppression = FindConnectedWeakEdges(non_max_suppression, strongEdgesRow(i),...
            strongEdgesCol(i));   
end
figure; imshow(non_max_suppression);
title('Edge Tracking'); 
            

% Find weak edges that are connected to strong edges 
function[non_max_suppression] = FindConnectedWeakEdges(non_max_suppression, row1, col1)
    for i = -3:1:3
        for j = -3:1:3
            if (row1+i > 0) && (col1+j > 0) && (row1+i < size(non_max_suppression,1)) && ...
                    (col1+j < size(non_max_suppression,2)) % Make sure we are not out of bounds
                if (non_max_suppression(row1+i,col1+j) > 0) && (non_max_suppression(row1+i,col1+j) < 1)
                    non_max_suppression(row1+i,col1+j) = 1;
                    non_max_suppression = FindConnectedWeakEdges(non_max_suppression, row1+i, col1+j);
                end
            end
        end
    end
end
    