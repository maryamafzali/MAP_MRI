function [RTOP,RTAP,RTPP,NG,NG_par,NG_per,PA_DTI, PA,delta_theta,I_trans,a_coeff,O,phi,kappa_iso,E, S0,D,V,FA,u0] = ...
    mapmri(S,bvals,bvecs,big_delta,little_delta,Nmax,method,n)

% big_delta = 39.1*10^-3;
% little_delta = 22.1*10^-3;
D0 = 0.003;

td = big_delta-little_delta/3;

q = zeros(3,size(bvecs,2));

M_sym = 1/6*(Nmax/2+1)*(Nmax/2+2)*(2*Nmax+3);
 
bvecs_1000 = bvecs(:,bvals<2000);
bvals_1000 = bvals(bvals<2000);
S1 = S(bvals<2000);

[Dcnls, ~] = constra_tensor_est(S1,bvals_1000,bvecs_1000);

[V, D] = eig(Dcnls);

D = 2*td*D;

[D, order] = sort(diag(D),'descend');  % sort eigenvalues in descending order
V = V(:,order);

ux = sqrt(D(1));
uy = sqrt(D(2));
uz = sqrt(D(3));

[ind , p_per , p0_per , a_par , a0_par , K , K_vert , K_par , B] = mapmri_p_k_a(Nmax);

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

zetax = ux/u0;
zetay = uy/u0;
zetaz = uz/u0;

for k = 1:size(bvecs,2)
    q(:,k) = V'*bvecs(:,k);
    q(:,k) = q(:,k).*repmat(sqrt(bvals(1,k)/td)/(2*pi),3,1);
end

phi_x = zeros(size(bvecs,2),size(ind,1));
phi_y = zeros(size(bvecs,2),size(ind,1));
phi_z = zeros(size(bvecs,2),size(ind,1));

for j = 1:size(ind,1)
    phi_x(:,j) = mapmri_phi(ind(j,1),ux,q(1,:));
    phi_y(:,j) = mapmri_phi(ind(j,2),uy,q(2,:));
    phi_z(:,j) = mapmri_phi(ind(j,3),uz,q(3,:));
end
phi = phi_x.*phi_y.*phi_z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constraining the estimation

switch method
    case 'MAP'
        grid_size = 20.4; %17;
        
        a = 1;
        %         for i = -grid_size:grid_size
        %             for j = -grid_size:grid_size
        %                 for k = 0:grid_size
        for i = -grid_size:0.8:grid_size
            for j = -grid_size:0.8:grid_size
                for k = 0:0.8:grid_size
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
        
        psi_x = zeros(size(samples,1),size(ind,1));
        psi_y = zeros(size(samples,1),size(ind,1));
        psi_z = zeros(size(samples,1),size(ind,1));
        
        w = -1*ones(size(samples,1),1);
        for i = 1:size(w)
            if samples(i,3) == 0
                w(i) = -0.5;
            end
        end
        
        for j = 1:size(ind,1)
            psi_x(:,j) = mapmri_psi(ind(j,1),sqrt(D(1)),samples(:,1));
            psi_y(:,j) = mapmri_psi(ind(j,2),sqrt(D(2)),samples(:,2));
            psi_z(:,j) = mapmri_psi(ind(j,3),sqrt(D(3)),samples(:,3));
        end
        
        psi = psi_x.*psi_y.*psi_z;
        
        S0 = mean(S(bvals<10));
        K1 = psi*delta_x^3;
        A = -[K1 ; w'*K1];
        b = [zeros(size(samples,1),1) ; 0.5];        
        
        phi0 = mean(phi(bvals<10,:),1);
        phi1 = phi(bvals>=10,:);
        phi1 = [phi0; phi1];
        H = 1*(phi1'*phi1);

        S1 = [S0 ; S(bvals>=10)];
        f = -1*phi1'*S1/S0;
        
        options = optimset('Algorithm','interior-point-convex','UseParallel',true);
        
        x = quadprog(sparse(H),f,A,b,[], [], [], [], [], options);
        
        a_coeff = x;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %MAPL
    case 'MAPL'
        
        U1 = zeros(3,size(ind,1),size(ind,1));
        T1 = zeros(3,size(ind,1),size(ind,1));
        S1 = zeros(3,size(ind,1),size(ind,1));
        U = zeros(size(ind,1));
        
        for i = 1:size(ind,1)
            for j = 1:size(ind,1)
                
                U1(1,i,j) = mapl_U_mex(ind(i,1),ind(j,1));
                U1(2,i,j) = mapl_U_mex(ind(i,2),ind(j,2));
                U1(3,i,j) = mapl_U_mex(ind(i,3),ind(j,3));
                
                T1(1,i,j) = mapl_T_mex(ind(i,1),ind(j,1));
                T1(2,i,j) = mapl_T_mex(ind(i,2),ind(j,2));
                T1(3,i,j) = mapl_T_mex(ind(i,3),ind(j,3));
                
                S1(1,i,j) = mapl_S_mex(ind(i,1),ind(j,1));
                S1(2,i,j) = mapl_S_mex(ind(i,2),ind(j,2));
                S1(3,i,j) = mapl_S_mex(ind(i,3),ind(j,3));
                
                U(i,j) = ux^3/(uy*uz)*S1(1,i,j)*U1(2,i,j)*U1(3,i,j) +2*ux*uy/uz*T1(1,i,j)*T1(2,i,j)*U1(3,i,j) + ...
                    uy^3/(ux*uz)*S1(2,i,j)*U1(3,i,j)*U1(1,i,j) + 2*uy*uz/ux*T1(2,i,j)*T1(3,i,j)*U1(1,i,j) + ...
                    uz^3/(ux*uy)*S1(3,i,j)*U1(1,i,j)*U1(2,i,j) + 2*ux*uz/uy*T1(1,i,j)*T1(3,i,j)*U1(2,i,j);
                
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
            S_lambda = phi*inv(phi'*phi+lrange(i)*U)*phi';
            gcvnew = norm(S - S_lambda*S)/(length(S)-trace(S_lambda));
        end
        lopt = lrange(i - 1);
        
        c = inv(phi'*phi + lopt*U)*phi'*S;
        S0 = mean(S(bvals<10));
        %         S0 = c'*B;
        a_coeff = c/S0;
        
    case 'FREE'
        c = inv(phi'*phi)*phi'*S;
        S0 = mean(S(bvals<10));
        a_coeff = c/S0;
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spherical part
ind_sph = zeros(round(M_sym),3);
p = 1;
for N = 0:2:Nmax
    for j = 1:Nmax
        for l = 0:2:Nmax
            for m = -l:l
                if l + 2*j - 2 == N
                    ind_sph(p,:) = [j l m];
                    p = p + 1;
                end
            end
        end
    end
end

E = zeros(size(bvecs,2),size(ind_sph,1));

for k = 1:size(ind_sph,1)
    E(:,k) = mapmri_Xi(ind_sph(k,1),ind_sph(k,2),ind_sph(k,3),u0,q,bvecs);
end

switch method
    case 'MAP'

        r = sqrt(samples(:,1).^2 + samples(:,2).^2 + samples(:,3).^2);
        r1 = r;
        r1(r1 == 0) = 1;
        svec = [samples(:,1)./r1 samples(:,2)./r1 samples(:,3)./r1];
        
        gamma = zeros(length(r),size(ind_sph,1));
        
        for j = 1:size(ind_sph,1)
            gamma(:,j) = (-1)^(ind_sph(j,1)-1)/(sqrt(2)*pi*u0^3)*(r.^2/2/u0^2).^(ind_sph(j,2)/2).*exp(-r.^2/2/u0^2)...
                .*mapmri_LaguerreL(ind_sph(j,1)-1,(ind_sph(j,2)+1/2),r.^2/u0^2).*spharm(ind_sph(j,2),ind_sph(j,3),svec);
        end
                
%         S0 = mean(S(bvals<10));
        K1 = gamma*delta_x^3;
        A = -[K1 ; w'*K1];
        b = [zeros(size(samples,1),1) ; 0.5];
        
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
        U1 = zeros(3,size(ind,1),size(ind,1));
        T1 = zeros(3,size(ind,1),size(ind,1));
        S1 = zeros(3,size(ind,1),size(ind,1));
        U = zeros(size(ind,1));
        
        for i = 1:size(ind,1)
            for j = 1:size(ind,1)
                
                U1(1,i,j) = mapl_U_mex(ind(i,1),ind(j,1));
                U1(2,i,j) = mapl_U_mex(ind(i,2),ind(j,2));
                U1(3,i,j) = mapl_U_mex(ind(i,3),ind(j,3));
                
                T1(1,i,j) = mapl_T_mex(ind(i,1),ind(j,1));
                T1(2,i,j) = mapl_T_mex(ind(i,2),ind(j,2));
                T1(3,i,j) = mapl_T_mex(ind(i,3),ind(j,3));
                
                S1(1,i,j) = mapl_S_mex(ind(i,1),ind(j,1));
                S1(2,i,j) = mapl_S_mex(ind(i,2),ind(j,2));
                S1(3,i,j) = mapl_S_mex(ind(i,3),ind(j,3));
                U(i,j) = ux^3/(uy*uz)*S1(1,i,j)*U1(2,i,j)*U1(3,i,j) +2*ux*uy/uz*T1(1,i,j)*T1(2,i,j)*U1(3,i,j) + ...
                    uy^3/(ux*uz)*S1(2,i,j)*U1(3,i,j)*U1(1,i,j) + 2*uy*uz/ux*T1(2,i,j)*T1(3,i,j)*U1(1,i,j) + ...
                    uz^3/(ux*uy)*S1(3,i,j)*U1(1,i,j)*U1(2,i,j) + 2*ux*uz/uy*T1(1,i,j)*T1(3,i,j)*U1(2,i,j);
                
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
        S0 = mean(S(bvals<10));
        kappa = c/S0;
        
    case 'FREE'
        
        c = inv(E'*E)*E'*S;
        S0 = mean(S(bvals<10));
        kappa = c/S0;
        
end

kappa_iso = zeros(size(ind,1),1);

for N = 0:2:Nmax
    kappa_iso((N/2+1)*(N/2+2)*(2*N+3)/6) = kappa((N/2+1)*(N/2+2)*(2*N+3)/6);
end

O = zeros(size(ind,1),1);

p = 1;
for N = 0:2:Nmax
    for i = 0:Nmax
        for j = 0:Nmax
            for k = 0:Nmax
                if i + j + k == N
                    if i/2 == round(i/2) && j/2 == round(j/2) && k/2 == round(k/2)
                        O(p) = kappa((N/2+1)*(N/2+2)*(2*N+3)/6)*...
                            sqrt(factorial(i)*factorial(j)*factorial(k))/...
                            (factd(i)*factd(j)*factd(k));
                    end
                    p = p + 1;
                end
            end
        end
    end
end

Tx = zeros(size(ind,1),size(ind,1));
Ty = zeros(size(ind,1),size(ind,1));
Tz = zeros(size(ind,1),size(ind,1));

for i = 1:size(ind,1)
    for j = 1:size(ind,1)
        if ind(i,1)/2 == round(ind(i,1)/2) && ind(j,1)/2 == round(ind(j,1)/2)
            Tx(i,j) = mapmri_T(ind(i,1),ind(j,1),zetax);
        end
        if ind(i,2)/2 == round(ind(i,2)/2) && ind(j,2)/2 == round(ind(j,2)/2)
            Ty(i,j) = mapmri_T(ind(i,2),ind(j,2),zetay);
        end
        if ind(i,3)/2 == round(ind(i,3)/2) && ind(j,3)/2 == round(ind(j,3)/2)
            Tz(i,j) = mapmri_T(ind(i,3),ind(j,3),zetaz);
        end
    end
end

T = Tx.*Ty.*Tz;

clear Tx Ty Tz

[X,Y,Z] = sphere(n);
xyz = [X(:),Y(:),Z(:)];
% sx = size(X);

I = zeros(size(xyz,1),1);

rho = [X(:)/ux, Y(:)/uy , Z(:)/uz]';

mag_rho = sqrt(rho(1,:).^2 + rho(2,:).^2 + rho(3,:).^2);

mag_rho1 = mag_rho;
mag_rho1(mag_rho1==0) = 1;

beta = 2*[rho(1,:)./mag_rho1; rho(2,:)./mag_rho1; rho(3,:)./mag_rho1];

for p = 1:length(mag_rho)
    I(p) = mapmri_ODF(ind,D,beta(:,p),mag_rho(p),a_coeff);
end

I_trans = [xyz(:,1).*I xyz(:,2).*I xyz(:,3).*I]*V';

RTOP = mapmri_RTOP(K,D,a_coeff);
RTAP = mapmri_RTAP(K_vert,sqrt(D(2)),sqrt(D(3)),a_coeff);
RTPP = mapmri_RTPP(K_par,sqrt(D(1)),a_coeff);
NG = mapmri_NG(a_coeff);
NG_par = mapmri_NG_par(a0_par,a_par,a_coeff);
NG_per = mapmri_NG_per(p0_per,p_per,a_coeff);

[PA , theta_PA] = mapmri_PA(O,T,D,u0,a_coeff);

[PA_DTI , theta_DTI] = mapmri_PA_DTI(ux,uy,uz,u0);

delta_theta = theta_PA - theta_DTI;

% MD = (D(1)+D(2)+D(3))/3;

R = Dcnls/trace(Dcnls);

FA = sqrt(0.5*(3-1/trace(R^2)));

end
