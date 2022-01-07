% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
function [rG, vG] = Gauss(Rho1, Rho2, Rho3, R1, R2, R3, t1, t2, t3, r2_true)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
% This function uses the Gauss method with iterative
% improvement (Algorithms 5.5 and 5.6) to calculate the state
% vector of an orbiting body from angles-only observations at
% three closely-spaced times.
%
% mu - the gravitational parameter (kmˆ3/sˆ2)
% t1, t2, t3 - the times of the observations (s)
% tau, tau1, tau3 - time intervals between observations (s)
% R1, R2, R3 - the observation site position vectors
% at t1, t2, t3 (km)
% Rho1, Rho2, Rho3 - the direction cosine vectors of the
% satellite at t1, t2, t3
% p1, p2, p3 - cross products among the three direction
% cosine vectors
% Do - scalar triple product of Rho1, Rho2 and
% Rho3
% D - Matrix of the nine scalar triple products
% of R1, R2 and R3 with p1, p2 and p3
% E - dot product of R2 and Rho2
% A, B - constants in the expression relating
% slant range to geocentric radius
% a,b,c - coefficients of the 8th order polynomial
% in the estimated geocentric radius x
% x - positive root of the 8th order polynomial
% rho1, rho2, rho3 - the slant ranges at t1, t2, t3
% r1, r2, r3 - the position vectors at t1, t2, t3 (km)
% r_old, v_old - the estimated state vector at the end of
% Algorithm 5.5 (km, km/s)
% rho1_old,
% rho2_old, and
% rho3_old - the values of the slant ranges at t1, t2,
% t3 at the beginning of iterative
% improvement (Algorithm 5.6) (km)
% diff1, diff2,
% and diff3 - the magnitudes of the differences between
% the old and new slant ranges at the end
% of each iteration
% tol - the error tolerance determining
% convergence
% n - number of passes through the
% iterative improvement loop
% nmax - limit on the number of iterations
% ro, vo - magnitude of the position and
% velocity vectors (km, km/s)
% vro - radial velocity component (km)
% a - reciprocal of the semimajor axis (1/km)
% v2 - computed velocity at time t2 (km/s)
% r, v - the state vector at the end of
% Algorithm 5.6 (km, km/s)
% ------------------------------------------------------------
global mu %katsayilar

%...Equations 5.98:
tau1 = t1 - t2;
tau3 = t3 - t2;
%...Equation 5.101:
tau = tau3 - tau1;
%...Independent cross products among the direction cosine vectors:
p1 = cross(Rho2,Rho3);
p2 = cross(Rho1,Rho3);
p3 = cross(Rho1,Rho2);
%...Equation 5.108:
Do = dot(Rho1,p1);
%...Equations 5.109b, 5.110b and 5.111b:
D = [[dot(R1,p1) dot(R1,p2) dot(R1,p3)]
    [dot(R2,p1) dot(R2,p2) dot(R2,p3)]
    [dot(R3,p1) dot(R3,p2) dot(R3,p3)]];
%...Equation 5.115b:
E = dot(R2,Rho2);
%...Equations 5.112b and 5.112c:
A = 1/Do*(-D(1,2)*tau3/tau + D(2,2) + D(3,2)*tau1/tau);
B = 1/6/Do*(D(1,2)*(tau3^2 - tau^2)*tau3/tau ...
    + D(3,2)*(tau^2 - tau1^2)*tau1/tau);
%...Equations 5.117:
a = -(A^2 + 2*A*E + norm(R2)^2);
b = -2*mu*B*(A + E);
c = -(mu*B)^2;

%...Calculate the roots of Equation 5.116 using MATLAB’s
% polynomial ‘roots’ solver:
Roots = roots([1 0 a 0 0 b 0 0 c]);

%...Find the positive real root:
    x = posroot(Roots, r2_true);

% katsayilar=[katsayilar; a b c x];    
%...Equations 5.99a and 5.99b:
f1 = 1 - 1/2*mu*tau1^2/x^3;
f3 = 1 - 1/2*mu*tau3^2/x^3;
%...Equations 5.100a and 5.100b:
g1 = tau1 - 1/6*mu*(tau1/x)^3;
g3 = tau3 - 1/6*mu*(tau3/x)^3;
%...Equation 5.112a:
rho2 = A + mu*B/x^3;
%...Equation 5.113:
rho1 = 1/Do*((6*(D(3,1)*tau1/tau3 + D(2,1)*tau/tau3)*x^3 ...
    + mu*D(3,1)*(tau^2 - tau1^2)*tau1/tau3) ...
    /(6*x^3 + mu*(tau^2 - tau3^2)) - D(1,1));
%...Equation 5.114:
rho3 = 1/Do*((6*(D(1,3)*tau3/tau1 - D(2,3)*tau/tau1)*x^3 ...
    + mu*D(1,3)*(tau^2 - tau3^2)*tau3/tau1) ...
    /(6*x^3 + mu*(tau^2 - tau3^2)) - D(3,3));
%...Equations 5.86:
r1 = R1 + rho1*Rho1;
r2 = R2 + rho2*Rho2;
r3 = R3 + rho3*Rho3;
%...Equation 5.118:
v2 = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1);
%...Save the initial estimates of r2 and v2:

rG = r2;
vG = v2;
return
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
% Subfunction used in the main body:
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
function x = posroot(Roots, r2_trueroot)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
%
% This subfunction extracts the positive real roots from
% those obtained in the call to MATLAB's 'roots' function.
% If there is more than one positive root, the user is
% prompted to select the one to use.
%
% x - the determined or selected positive root
% Roots - the vector of roots of a polynomial
% posroots - vector of positive roots
%
% User M-functions required: none
% ------------------------------------------------------------
%...Construct the vector of positive real roots:
posroots = Roots(find(Roots>0 & ~imag(Roots)));
npositive = length(posroots);
%...Exit if no positive roots exist:
if npositive == 0
    fprintf('\n\n ** There are no positive roots. \n\n')
    return
end
%...If there is more than one positive root, output the
%...roots to the command window and prompt the user to
%...select which one to use:
if npositive == 1
    x = posroots(1);
else
   x = posroots(1);
  for jj=2:npositive;   
      if abs(posroots(jj)-r2_trueroot)<abs(x-r2_trueroot)
          x = posroots(jj);
      end
  end
end
return
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜