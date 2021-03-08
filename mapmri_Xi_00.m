function Xi = mapmri_Xi_00(J,u0,mag_q)

% mag_q = sqrt(q(1,:).^2 + q(2,:).^2 + q(3,:).^2);
% mag_q = mag_q';

% mag_q1 = mag_q;
% mag_q1(mag_q1 == 0) = 1;
% bvecs = [q(1,:)'./mag_q1 q(2,:)'./mag_q1 q(3,:)'./mag_q1];

Xi = sqrt(4*pi).*exp(-2*pi^2*u0^2*mag_q.^2).*mapmri_LaguerreL(J-1,1/2,4*pi^2*u0^2*mag_q.^2)*1/sqrt(4*pi);

end
