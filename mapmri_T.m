function T = mapmri_T(m,n,zeta)

coeff = 0;

for r = 0:2:m
    for s = 0:2:n
        coeff = coeff + zeta^(n-s)*((1+zeta^2)/2)^-((m+n-r-s+1)/2)*(-1)^((r+s)/2)...
            *factd(m+n-r-s-1)/factorial(m-r)/factorial(n-s)/factd(r)/factd(s);
    end
end

T = sqrt(factorial(m)*factorial(n))*coeff;

end