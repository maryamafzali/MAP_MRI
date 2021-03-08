function [ind , p_per , p0_per , a_par , a0_par , K , K_vert , K_par, B] = mapmri_p_k_a(Nmax)

M_sym = 1/6*(Nmax/2+1)*(Nmax/2+2)*(2*Nmax+3);

ind = zeros(round(M_sym),3);

K = zeros(round(M_sym),1);
B = zeros(round(M_sym),1);

K_vert = zeros(round(M_sym),1);
K_par = zeros(round(M_sym),1);

p_per = zeros(round(M_sym),1);
p0_per = zeros(round(M_sym),1);

a_par = zeros(round(M_sym),1);
a0_par = zeros(round(M_sym),1);

p = 1;
for N = 0:2:Nmax
    for k = 0:Nmax
        for j = 0:Nmax
            for i = 0:Nmax
                if i + j + k == N
                    ind(p,:) = [i j k];
                    if i/2 == round(i/2)
                        p_per(p) = (-1)^(i/2)*sqrt(factorial(i))/(factd(i));
                        if j == 0 && k == 0
                            p0_per(p) = (-1)^(i/2)*sqrt(factorial(i))/(factd(i));
                        end
                    end
                    if j/2 == round(j/2) && k/2 == round(k/2)
                        a_par(p) = (-1)^((j+k)/2)*sqrt(factorial(j)*factorial(k))/(factd(j)*factd(k));
                        if i == 0
                            a0_par(p) = (-1)^((j+k)/2)*sqrt(factorial(j)*factorial(k))/(factd(j)*factd(k));
                        end
                    end
                    if i/2 == round(i/2) && j/2 == round(j/2) && k/2 == round(k/2)
                        K(p) = (-1)^((i+j+k)/2)*sqrt(factorial(i)*factorial(j)*factorial(k))...
                            /(factd(i)*factd(j)*factd(k));
                        B(p) = sqrt(factorial(i)*factorial(j)*factorial(k))...
                            /(factd(i)*factd(j)*factd(k));
                        K_vert(p) = (-1)^((j+k)/2)*sqrt(factorial(i)*factorial(j)*factorial(k))...
                            /(factd(i)*factd(j)*factd(k));
                        K_par(p) = (-1)^(i/2)*sqrt(factorial(i)*factorial(j)*factorial(k))/(factd(i)*factd(j)*factd(k));
                    end
                    p = p + 1;
                end
            end
        end
    end
end