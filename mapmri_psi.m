function psi = mapmri_psi(n,u,x)

psi = exp(-x.^2/(2*u^2)).*hermite(n,x/u)/(sqrt(2^(n+1)*pi*factorial(n))*u);

end