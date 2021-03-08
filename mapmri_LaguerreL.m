function result = mapmri_LaguerreL(n, alpha, x)

if n == 0
    result = 1;
else
    Hn1 = 1;
    result = -x + (alpha + 1);
    
    for i = 2:n
        Hn2=Hn1;
        Hn1=result;
        result=(((2.0 * i + alpha - 1.0)-x).*Hn1 - ( i + alpha - 1.0) .* Hn2) / i;
    end
end

