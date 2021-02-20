function [sobXxX] = custom_sobel(kernel_size)
% kernel_size must be odd
% from https://stackoverflow.com/a/57399675/9275698
halfway = floor(kernel_size/2);
step = -halfway : halfway;

i_matrix = repmat(step,[kernel_size 1]);
j_matrix = i_matrix';

sobXxX = i_matrix ./ ( i_matrix.*i_matrix + j_matrix.*j_matrix );
sobXxX(halfway+1,halfway+1) = 0; % deals with NaN in middle