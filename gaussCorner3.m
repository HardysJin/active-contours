function [result] = gaussCorner3(I, corners, cornerSize, minSigma, maxSigma, steps)
Y = cornerSize(1);
X = cornerSize(2);
% [Y, X] = size(I);
ratio = Y/X;
% tmp = zeros(100, 100);
% tmp(50:100, 50:100) = 1;
% imshow(tmp)
sigmaList = minSigma:steps:maxSigma;
s = size(sigmaList, 2);
result = zeros(size(I));
roi = zeros([size(I), s]);
% mid = floor(size(I) / 2);

interval = 3/s; % Standard Normal Distribution cdf step

for i = 1:s
    sigma = sigmaList(i);
    
    pro = normcdf([-3, -3+interval*i]);
    
    p = (pro(2) - pro(1)) * 2;
    
    % x*y = area; y/x = ratio => y^2 = area*ratio
    area = Y * X * p;
    
    roi_y = sqrt(area * ratio);
    roi_x = area / roi_y;
    
    x_t = floor(roi_x/2);
    y_t = floor(roi_y/2);
    
    z = zeros(size(I));
    for j = 1:size(corners, 1)
        corner = corners(j,:);
        x_l = max(corner(1) - x_t, 1);
        x_r = min(corner(1) + x_t, size(I,2));
        
        y_l = max(corner(2) - y_t, 1);
        y_r = min(corner(2) + y_t, size(I,1));
        
        z(y_l:y_r, x_l:x_r) = 1;
    end

    roi(:, :, i) = z;

    if i > 1
        z = z - roi(:,:,i-1);
    end
    % imshow(z);
    % https://www.mathworks.com/matlabcentral/answers/154064-what-do-fspecial-image-filtering-parameters-mean
    hsize = 2*ceil(2.6*sigma)+1;
    h = fspecial('gaussian', hsize, sigma);
    J = roifilt2(h,I,z);

    if size(J,1) > 0
        result = result + z .* J;
    end
%     if i == 2
%         get3D(roi .* imgaussfilt(tmp, sigma));
%         error("2");
%     end
%     get3D(result);
%     error("1");
end

assignin('base','roi',roi(:,:,s));

z = ones(size(I));
z = z - roi(:,:,s);

hsize = 2*ceil(2.6*maxSigma)+1;
h = fspecial('gaussian', hsize, maxSigma);
J = roifilt2(h,I,z);

result = result + z .* J;

% get3D(result);