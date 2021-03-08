function T = mapl_T(n,m)

T = (-1)^(n+1)*pi^(1.5)*(kronDel(m-n)*(1+2*n)+ ...
    kronDel(m+2-n)*sqrt(n*(n-1)) + ...
    kronDel(n+2-m)*sqrt(m*(m-1)));

end

