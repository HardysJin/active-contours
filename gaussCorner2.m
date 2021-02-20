function [result] = gaussCorner2(I, minSigma, maxSigma, steps)

[Y, X] = size(I);
ratio = Y/X;
% tmp = zeros(100, 100);
% tmp(50:100, 50:100) = 1;
% imshow(tmp)
sigmaList = minSigma:steps:maxSigma;
s = size(sigmaList, 2);
result = zeros(Y, X);
roi = zeros([Y, X, s]);
mid = floor(size(I) / 2);

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
    
    x_l = max(mid(2) - x_t, 1);
    x_r = min(mid(2) + x_t, size(I,2));
    
    y_l = max(mid(1) - y_t, 1);
    y_r = min(mid(1) + y_t, size(I,1));
    
    z = zeros(size(I));
    z(y_l:y_r, x_l:x_r) = 1;
    imshow(z)
    roi(:, :, i) = z;
    
    if i > 1
        z = z - roi(:,:,i-1);
    end
    % get3D(z);
    
    result = result + z .* imgaussfilt(I, sigma);
%     if i == 2
%         get3D(roi .* imgaussfilt(tmp, sigma));
%         error("2");
%     end
%     get3D(result);
%     error("1");
end
assignin('base','Eedge',roi(:,:,s));

z = ones(size(I));
z = z - roi(:,:,s);
result = result + z .* imgaussfilt(I, maxSigma);

% get3D(result);