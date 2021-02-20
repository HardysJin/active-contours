function [Ainv] = getInternalEnergyMatrixBonus(nPoints, alpha, beta, gamma)

Xi2 = beta;
Xi1 = -alpha - 4*beta;
Xi  = 2*alpha + 6*beta;

A1 = Xi * eye(nPoints);

A2 = Xi1 * eye(nPoints);
% circshift(A2, -1)
A2 = [A2(:,end) A2(:,1:end-1)]; % shift right one;

A3 = Xi1 * eye(nPoints);
A3 = [A3(:,2:end) A3(:,1)]; % shift left one;

A4 = Xi2 * eye(nPoints);
A4 = [A4(:,end-1:end) A4(:,1:end-2)]; % shift right two;

A5 = Xi2 * eye(nPoints);
A5 = [A5(:,3:end) A5(:,1:2)]; % shift left two;

A = A1 + A2 + A3 + A4 + A5;

% assignin('base', 'A', A);
Ainv=inv(A + gamma * eye(nPoints)); 

end