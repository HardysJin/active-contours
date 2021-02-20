clear;
close all;

% Parameters (play around with different images and different parameters)
N=600;
Wedge=4;
Wline=1;
Wterm=0.1;

Wconcave=1;
Wconvex=1.5;

kappa=0.05;
sigma=10;
% minSigma=7.6;
minSigma=7.6;
kernel = 9;

alpha=0.2;
beta=0.1;
gamma=0.1;

% Load image
img = 'shape';
suff = '.png';
I = imread("images/" + img + suff);
I = im2double(I);
if (ndims(I) == 3)
    I = rgb2gray(I);
end

switch img
    case {'circle', 'square'}
        I = imcomplement(I);
end
    

% Initialize the snake
[x, y] = initializeSnake(I);

% Calculate external energy
I = double(I);
I_smooth = imgaussfilt(I, sigma);

% --------------------------------------------------------------------
% Corners definitions:
corners = [  185.8750 76.3750; 
            99.8750 76.375; 
            156.3750 185.3750;
            97.3750 271.3750;
            228.8750 271.3750];

if strcmp(img,'shape')
    sigma=10;
    kappa=0.05;
    Wterm=0.1;
    Wconcave=11;
    Wconvex=6;
    cornerSize = [61 61];
    result = gaussCorner3(I, corners, cornerSize, minSigma, sigma, 0.1);
    I_smooth = result;
end

% get3D(result);
% get3D(imgradient(result));
% error("1")
% --------------------------------------------------------------------


Eext = getExternalEnergy(I_smooth,Wline,Wedge,Wterm,Wconcave,Wconvex,kernel);

% Calculate matrix A^-1 for the iteration
Ainv = getInternalEnergyMatrixBonus(size(x,2), alpha, beta, gamma);

% Get fx and fy - known gradient Eext

% [fx, fy] = imgradientxy(Eext, 'sobel');
% 
% A = [ 1 2 1 ]' * [1 2 1];
% sob3x3 = [ 1 2 1 ]' * [1 0 -1];
% sob5x5 = conv2( A, sob3x3 );
% sob7x7 = conv2( A, sob5x5 );
% sob7x7 = conv2( A, sob7x7 );
% sob7x7 = sob7x7/max(sob7x7(:));
% sob7x7  = [ 3/18 2/13 1/10 0  -1/10 -2/13 -3/18;
%             3/13 2/8  1/5  0  -1/5  -2/8  -3/13;
%             3/10 2/5  1/2  0  -1/2  -2/5  -3/10;
%             3/9  2/4  1/1  0  -1/1  -2/4  -3/9;
%             3/10 2/5  1/2  0  -1/2  -2/5  -3/10;
%             3/13 2/8  1/5  0  -1/5  -2/8  -3/13;
%             3/18 2/13 1/10 0  -1/10 -2/13 -3/18;];
%         
% sob9x9 = [  4/32 3/25 2/20 1/17  0/16  -1/17  -2/20  -3/25  -4/32;
%             4/25 3/18 2/13 1/10   0/9  -1/10  -2/13  -3/18  -4/25;
%             4/20 3/13  2/8  1/5   0/4   -1/5   -2/8  -3/13  -4/20;
%             4/17 3/10  2/5  1/2   0/1   -1/2   -2/5  -3/10  -4/17;
%             4/16  3/9  2/4  1/1     0   -1/1   -2/4   -3/9  -4/16;
%             4/17 3/10  2/5  1/2   0/1   -1/2   -2/5  -3/10  -4/17;
%             4/20 3/13  2/8  1/5   0/4   -1/5   -2/8  -3/13  -4/20;
%             4/25 3/18 2/13 1/10   0/9  -1/10  -2/13  -3/18  -4/25;
%             4/32 3/25 2/20 1/17  0/16  -1/17  -2/20  -3/25  -4/32];
sobelX = -custom_sobel(kernel);
fx = conv2(Eext, sobelX, 'same');
fy = conv2(Eext, sobelX', 'same');


imshow(I, 'InitialMagnification', 300); hold on;
[xx,yy]=ndgrid(1:3:size(fx,2),1:3:size(fy,1));
quiver(xx,yy, -10*kappa*fx(1:3:end,1:3:end)', -10*kappa*fy(1:3:end,1:3:end)');
title('The external force')
% hold on;

x = x';
y = y';
% Iterate and update positions
displaySteps = floor(N/10);
figure;
for i=1:N
    % Iterate
    [x,y] = iterate(Ainv, x, y, Eext, gamma, kappa, fx, fy);
    % size(x)
    % Plot intermediate result
    imshow(I,'InitialMagnification',300); 
    hold on;
    % plot([x; x(1)], [y; y(1)], 'r');
    plot([x; x(1)], [y; y(1)], 'r--o'); 

    % Display step
    if(mod(i,displaySteps)==0)
        fprintf('%d/%d iterations\n',i,N);
    end
    
    pause(0.0000001)
end
 
if(displaySteps ~= N)
    fprintf('%d/%d iterations\n',N,N);
end