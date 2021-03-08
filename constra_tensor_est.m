function [Dcnls, S0, resnorm] = constra_tensor_est(S,bvals,bvecs)
% S is a mx1 vector
% bvals is 1xm
m = size(S,1);

S = double(S);
weight = S;
weight(weight<=0) = eps;
y = log(weight);

% y = log(S);

b = repmat(bvals',1,6);

g = bvecs';

G = [-g(:,1).^2, -g(:,2).^2, -g(:,3).^2, -2*g(:,1).*g(:,2), -2*g(:,2).*g(:,3), -2*g(:,1).*g(:,3)];

W = G.*b;

W = [ones(m,1),W];

ungama = lscov(W,y);

gama = lscov(W,y,weight);

A = [gama(2),gama(5),gama(7);gama(5),gama(3),gama(6);...
    gama(7),gama(6),gama(4)];
Dwlls = A;

Dlls = [ungama(2),ungama(5),ungama(7);ungama(5),ungama(3),ungama(6);...
    ungama(7),ungama(6),ungama(4)];

U0 = cholmod(A);

u0 = [log(S(1)) U0(1,1) U0(1,2) U0(1,3) U0(2,2) U0(2,3) U0(3,3)]';

options = optimset('TolFun',0.001,...
                   'Display','off',...
                   'Algorithm','levenberg-marquardt');
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@difsig,u0,W,S,[],[],options);

A1 = [x(2) x(3) x(4); 0 x(5) x(6); 0 0 x(7)];
Dcnls = A1'*A1;
S0 = exp(x(1));