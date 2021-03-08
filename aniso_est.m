function [trace_L, FA, cl, cp, cs, cm,D] = aniso_est(Dcnls)

[V, D] = eig(Dcnls);
[D, order] = sort(diag(D),'descend');  %# sort eigenvalues in descending order
V = V(:,order);

FA = sqrt((D(1) - D(2))^2 + (D(2) - D(3))^2 + (D(1) - D(3))^2)/sqrt(D(1)^2 + D(2)^2 + D(3)^2)/sqrt(2);

trace_L = trace(D(1) + D(2) + D(3));

cl = (D(1) - D(2)) / D(1);

cp = (D(2) - D(3)) / D(1);

cs = D(3) / D(1);

V = abs(V);       %get the x,y,z components of e3
r = V(1,1);
g = V(2,1);
b = V(3,1);

rgb = cat(3,r,g,b);
FAmap3 = cat(3,FA,FA,FA); %needed to combine with the colormap, need 3D.
cm = rgb.*FAmap3; %FA weighting for colormap
% figure, image (cm);axis image





