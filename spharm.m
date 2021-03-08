% This function generates the Spherical Harmonics basis functions of degree
% L and order M.
%
% SYNTAX: [Ymn,THETA,PHI,X,Y,Z]=spharm4(L,M,RES,PLOT_FLAG);
%
% INPUTS:
%
% L         - Spherical harmonic degree, [1x1]
% M         - Spherical harmonic order,  [1x1]
% RES       - Vector of # of points to use [#Theta x #Phi points],[1x2] or [2x1]
% PLOT_FLAG - Binary flag to generates a figure of the spherical harmonic surfaces (DEFAULT=1)
%
%
% OUTPUTS:
%
% Ymn   - Spherical harmonics coordinates, [RES(1) x RES(2)]
% THETA - Circumferential coordinates, [RES(1) x RES(2)]
% PHI   - Latitudinal coordinates, [RES(1) x RES(2)]
% X,Y,Z - Cartesian coordinates of magnitude, squared, spherical harmonic surface points, [RES(1) x RES(2)]
%
%
% NOTE: It is very important to keep the various usages of THETA and PHI
% straight.  For this function THETA is the Azimuthal/Longitude/Circumferential
% coordinate and is defined on the interval [0,2*pi], whereas PHI is the
% Altitude/Latitude/Elevation and is defined on the interval [0,pi].  Also note that
% the conversion to cartesian coordinates requires that PHI be offset by pi/2 so
% that the conversion is on the interval [-pi/2,pi/2].
%
% DBE 2005/09/30

function [re_Ymn]=spharm(L,M,bvec)

if L<M, error('The ORDER (M) must be less than or eqaul to the DEGREE(L).'); end

anglesOut = cart2sph_bis(bvec);
THETA = anglesOut(:,1);
PHI = anglesOut(:,2);

Lmn=legendre(L,cos(THETA));

if L~=0
    if M>0
        Lmn=Lmn(M+1,:);
    end
    if M == 0
        Lmn=Lmn(M+1,:);
    end
    %     if M<0
    %         Lmn=Lmn(-M+1,:)*(-1)^M*factorial(L+M)/factorial(L-M);
    %     end
end

a1=((2*L+1)/(4*pi));
a2=factorial(L-M)/factorial(L+M);
C=sqrt(a1*a2);

Ymn=C*Lmn(1,:)'.*exp(1i*M*PHI);

if M>0
    re_Ymn = sqrt(2)*(-1)^(M+1)*imag(Ymn);
end

if M == 0
    re_Ymn = Ymn;
end

if M<0
    a2=factorial(L+M)/factorial(L-M);
    C=sqrt(a1*a2);
    Lmn = Lmn(-M+1,:);
    re_Ymn=sqrt(2)*real(C*Lmn'.*exp(1i*-M*PHI));
end
return