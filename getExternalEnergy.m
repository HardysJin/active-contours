function [Eext] = getExternalEnergy(I,Wline,Wedge,Wterm,Wconscave,Wconvex,kernel)

% Eline = I(x,y)
Eline = I;

% Eedge = −|∇I(x,y)|^2 where |∇I(x,y)|^2 = |Ix|^2 + |Iy|^2

% [Gmag,Gdir] = imgradient(I, 'sobel');
% Eedge = -Gmag;


% A = [ 1 2 1 ]' * [1 2 1];
% sob3x3 = [ 1 2 1 ]' * [1 0 -1];
% sob5x5 = conv2( A, sob3x3 );
% sob7x7 = conv2( A, sob5x5 );
% % sob9x9 = conv2( A, sob7x7 );
% sob7x7 = sob7x7/max(sob7x7(:));
% sob7x7  = [ 3/18 2/13 1/10 0  -1/10 -2/13 -3/18;
%             3/13 2/8  1/5  0  -1/5  -2/8  -3/13;
%             3/10 2/5  1/2  0  -1/2  -2/5  -3/10;
%             3/9  2/4  1/1  0  -1/1  -2/4  -3/9;
%             3/10 2/5  1/2  0  -1/2  -2/5  -3/10;
%             3/13 2/8  1/5  0  -1/5  -2/8  -3/13;
%             3/18 2/13 1/10 0  -1/10 -2/13 -3/18;];
% % Ix = conv2(I, sob7x7, 'same');
% % Iy = conv2(I, sob7x7', 'same');
% % Eedge = - sqrt(Ix.^2 + Iy.^2);
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

sobolX = -custom_sobel(kernel);
Ix = conv2(I, sobolX, 'same');
Iy = conv2(I, sobolX', 'same');
Eedge = - sqrt(Ix.^2 + Iy.^2);

Ixx = conv2(Ix, sobolX, 'same');
Ixy = conv2(Ix, sobolX', 'same');
Iyy = conv2(Iy, sobolX', 'same');


% Eterm

% [Ix,Iy] = imgradientxy(I, 'sobel');
% [Ixx,Ixy] = imgradientxy(Ix, 'sobel');
% [~,Iyy] = imgradientxy(Iy, 'sobel');

% combine
Eterm = (Iyy.*Ix.^2 -2*Ixy.*Ix.*Iy + Ixx.*Iy.^2)./((1+Ix.^2 + Iy.^2).^(3/2));

Eterm(Eterm > 0) = Eterm(Eterm > 0) * Wconscave;
Eterm(Eterm < 0) = Eterm(Eterm < 0) * Wconvex;

% assignin('base','Eline',Eline);
% assignin('base','Eedge',Eedge);
% assignin('base','Eterm',Eterm);

% Eext
Eext = Wline*Eline + Wedge*Eedge - Wterm*abs(Eterm);

end

