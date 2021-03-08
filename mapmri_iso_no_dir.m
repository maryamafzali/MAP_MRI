function [kappa,E, u0] = ...
    mapmri_iso_no_dir(S,bvals,big_delta,little_delta,Nmax,method)

% big_delta = 39.1*10^-3;
% little_delta = 22.1*10^-3;
D0 = 0.003;

td = big_delta-little_delta/3;

% q = zeros(3,size(bvecs,2));

% M_sym = 1/6*(Nmax/2+1)*(Nmax/2+2)*(2*Nmax+3);
% ind = zeros(round(M_sym),3);

% p = 1;
% for N = 0:2:Nmax
%     for i = 0:Nmax
%         for j = 0:Nmax
%             for k = 0:Nmax
%                 if i + j + k == N
%                     ind(p,:) = [i j k];
%                     p = p + 1;
%                 end
%             end
%         end
%     end
% end

ind = [0,0,0; 2, 0, 0; 4, 0, 0; 6, 0, 0; 8, 0, 0];

% bvecs_1000 = bvecs;
% bvals_1000 = bvals;
% S1 = S;

% bvecs_1000 = bvecs(:,bvals<2000);
% bvals_1000 = bvals(bvals<2000);
% S1 = S(bvals<2000);

% [Dcnls, ~] = constra_tensor_est(S1,bvals_1000,bvecs_1000);
[D, S0, resnorm] = est_D(S,bvals);

% [V, D] = eig(Dcnls);

D = 2*td*D;

% [D, order] = sort(diag(D),'descend');  % sort eigenvalues in descending order
% V = V(:,order);
% 
u0 = sqrt(D);
% uy = sqrt(D(2));
% uz = sqrt(D(3));
% 
% a = -3;
% X = ux^2;
% Y = uy^2;
% Z = uz^2;
% 
% b = -(X + Y +Z);
% c = X*Y + X*Z + Y*Z;
% d = 3*X*Y*Z;
% 
% e = -b^3/27/a^3 + b*c/6/a^2 -d/2/a;
% f = c/3/a - b^2/9/a^2;
% 
% U = (e + sqrt(e^2 + f^3))^(1/3) + (e - sqrt(e^2 + f^3))^(1/3) -b/(3*a);
% 
% u0 = sqrt(U);

for k = 1:length(bvals)
    mag_q(k) = sqrt(bvals(1,k)/td)/(2*pi);
end

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

E = zeros(length(mag_q),Nmax/2 + 1);

for k = 1:Nmax/2 + 1
    E(:,k) = mapmri_Xi_00(k,u0,mag_q);
end

switch method
    case 'MAP'
        grid_size = 17;
        
        a = 1;
        for i = -grid_size:4:grid_size
            for j = -grid_size:4:grid_size
                for k = 0:4:grid_size
                    if norm([i j k])<=grid_size
                        samples(a,:) = [i j k];
                        a = a + 1;
                    end
                end
            end
        end
        rmax = sqrt(10*D0*td);
        delta_x = rmax/grid_size;
        samples = samples*delta_x;
        r = sqrt(samples(:,1).^2 + samples(:,2).^2 + samples(:,3).^2);
        r1 = r;
        r1(r1 == 0) = 1;
        svec = [samples(:,1)./r1 samples(:,2)./r1 samples(:,3)./r1];
        
        gamma = zeros(length(r),Nmax/2 + 1);
        
        for j = 1:Nmax/2 + 1
            gamma(:,j) = (-1)^(j-1)/(sqrt(2)*pi*u0^3).*exp(-r.^2/2/u0^2)...
                .*laguerreL(j-1,1/2,r.^2/u0^2)*0.2821;
        end
        
        w = -1*ones(size(samples,1),1);
        for i = 1:size(w)
            if samples(i,3) == 0
                w(i) = -0.5;
            end
        end
        
         S0 = mean(S(bvals<10));
% S0 = 1;
        K1 = gamma*delta_x^3;
        A = -[K1 ; w'*K1];
        b = [zeros(size(samples,1),1) ; 0.5];
        
%         c = inv(phi'*phi)*phi'*S;
%         S0 = c'*B;
        
%         E1 = E;
%         H = 1*(E1'*E1);
%         f = -1*E1'*S/S0;
        E0 = mean(E(bvals<10,:),1);
        E1 = E(bvals>=10,:);
        E1 = [E0; E1];
        H = 1*(E1'*E1);
        S1 = [S0 ; S(bvals>=10)];
        f = -1*E1'*S1/S0;
        
        options = optimset('Algorithm','interior-point-convex','UseParallel',true);
        
        x = quadprog(sparse(H),f,A,b,[], [], [], [], [], options);
        
        kappa = x;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %MAPL
    case 'MAPL'
        U1 = zeros(3,Nmax/2+1,Nmax/2+1);
        T1 = zeros(3,Nmax/2+1,Nmax/2+1);
        S1 = zeros(3,Nmax/2+1,Nmax/2+1);
        U = zeros(Nmax/2+1);
        
        for i = 1:Nmax/2+1
            for j = 1:Nmax/2+1
                
                U1(1,i,j) = mapl_U(2*i-2,2*j-2);
                U1(2,i,j) = 1/sqrt(4*pi);
                U1(3,i,j) = 1/sqrt(4*pi);
                
                T1(1,i,j) = mapl_T_mex(2*i-2,2*j-2);
                T1(2,i,j) = -5.5683;
                T1(3,i,j) = -5.5683;
                
                S1(1,i,j) = mapl_S_mex(2*i-2,2*j-2);
                S1(2,i,j) = -5.5683;
                S1(3,i,j) = -5.5683;
                U(i,j) = u0*S1(1,i,j)*U1(2,i,j)*U1(3,i,j) +2*u0*T1(1,i,j)*T1(2,i,j)*U1(3,i,j) + ...
                    u0*S1(2,i,j)*U1(3,i,j)*U1(1,i,j) + 2*u0*T1(2,i,j)*T1(3,i,j)*U1(1,i,j) + ...
                    u0*S1(3,i,j)*U1(1,i,j)*U1(2,i,j) + 2*u0*T1(1,i,j)*T1(3,i,j)*U1(2,i,j);
                
            end
        end
        lrange = 0.05:1/20:1;
        samp = size(lrange',1);
        gcvold = 10e10;
        gcvnew = gcvold;
        i = 1;
        while gcvold >= gcvnew && i < samp - 1
            gcvold = gcvnew;
            i = i + 1;
            S_lambda = E*inv(E'*E+lrange(i)*U)*E';
            gcvnew = norm(S - S_lambda*S)/(length(S)-trace(S_lambda));
        end
        lopt = lrange(i - 1);
        
        c = inv(E'*E + lopt*U)*E'*S;
%         S0 = c'*B;
         S0 = mean(S(bvals<10));
% S0 = 1;
        kappa = c/S0;
        
    case 'FREE'
        
        c = inv(E'*E)*E'*S;
%         S0 = c'*B;
         S0 = mean(S(bvals<10));
% S0 = 1;
        kappa = c/S0;
        
end

% kappa_iso = zeros(size(ind,1),1);
% 
% for N = 0:2:Nmax
%     kappa_iso((N/2+1)*(N/2+2)*(2*N+3)/6) = kappa((N/2+1)*(N/2+2)*(2*N+3)/6);
% end
end
