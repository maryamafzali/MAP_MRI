function phi = mapmri_phi(n,u,q)
phi = sqrt(-1)^-n/sqrt(2^n*factorial(n))*exp(-2*pi^2*q.^2*u^2).*hermite(n,2*pi*u*q);
end