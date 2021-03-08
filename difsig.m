function F = difsig(u,W)

% ro = [log(S(1)); U0(1,1)^2; U0(1,2)^2 + U0(2,2)^2; ...
%     U0(1,3)^2 + U0(2,3)^2 + U0(3,3)^2; U0(1,1)*U0(1,2); ...
%     U0(1,2)*U0(1,3) + U0(2,2)*U0(2,3); U0(1,1)*U0(1,3)];

% u0 = [log(S(1)) U0(1,1) U0(1,2) U0(1,3) U0(2,2) U0(2,3) U0(3,3)]';

U = [u(1); u(2)^2; u(3)^2 + u(5)^2; ...
    u(4)^2 + u(6)^2 + u(7)^2; u(2)*u(3); ...
    u(3)*u(4) + u(5)*u(6); u(2)*u(4)];

F = exp(W*U);