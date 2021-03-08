function Xi = mapmri_Xi(J,L,M,u0,q,bvecs)

mag_q = sqrt(q(1,:).^2 + q(2,:).^2 + q(3,:).^2);
mag_q = mag_q';
        
Xi = sqrt(4*pi)*1i^(-L)*(2*pi^2*u0^2*mag_q.^2).^(L/2).*exp(-2*pi^2*u0^2*mag_q.^2).*...
                    mapmri_LaguerreL(J-1,(L+1/2),4*pi^2*u0^2*mag_q.^2).*spharm(L,M,bvecs');

end



