function [newX, newY, fx, fy] = iterate(Ainv, x, y, Eext, gamma, kappa, fx, fy)

% x=min(max(x,1),size(fx,1));
% y=min(max(y,1),size(fx,2));

% Iterate
Fx = kappa*interp2(fx, x, y);
Fy = kappa*interp2(fy, x, y);
Fx(isnan(Fx))=0;
Fy(isnan(Fy))=0;

% more force to make snake rotate:
% assignin('base','Fx',Fx);
% assignin('base','Fy',Fy);



newX = Ainv * (gamma*x - Fx);
newY = Ainv * (gamma*y - Fy);





% Clamp to image size
[h, w] = size(Eext); % image height & width
newX = max(newX, 1);
newX = min(newX, w);
newY = max(newY, 1);
newY = min(newY, h);

% newX=min(max(newX,1),size(fx, 1));
% newY=min(max(newY,1),size(fx, 2));