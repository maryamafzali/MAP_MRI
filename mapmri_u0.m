function u0 = mapmri_u0(S,bvals,bvecs,td)

% big_delta = 39.1*10^-3;
% little_delta = 22.1*10^-3;
% D0 = 0.003;

bvecs_1000 = bvecs;
bvals_1000 = bvals;
S1 = S;

% bvecs_1000 = bvecs(:,bvals<2000);
% bvals_1000 = bvals(bvals<2000);
% S1 = S(bvals<2000);
 
[Dcnls, ~] = constra_tensor_est(S1,bvals_1000,bvecs_1000);

[~, D] = eig(Dcnls);

D = 2*td*D;

[D] = sort(diag(D),'descend');  % sort eigenvalues in descending order
% V = V(:,order);

ux = sqrt(D(1));
uy = sqrt(D(2));
uz = sqrt(D(3));

a = -3;
X = ux^2;
Y = uy^2;
Z = uz^2;

b = -(X + Y +Z);
c = X*Y + X*Z + Y*Z;
d = 3*X*Y*Z;

e = -b^3/27/a^3 + b*c/6/a^2 -d/2/a;
f = c/3/a - b^2/9/a^2;

U = (e + sqrt(e^2 + f^3))^(1/3) + (e - sqrt(e^2 + f^3))^(1/3) -b/(3*a);

u0 = sqrt(U);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spherical part
% ind_sph = zeros(round(M_sym),3);
% p = 1;
% for N = 0:2:Nmax
%     for j = 1:Nmax
%         for l = 0:2:Nmax
%             for m = -l:l
%                 if l + 2*j - 2 == N
%                     ind_sph(p,:) = [j l m];
%                     p = p + 1;
%                 end
%             end
%         end
%     end
% end
