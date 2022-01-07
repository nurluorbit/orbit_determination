function [Robject2, Vobject2] = Laplace(rhoHat1, rhoHat2, rhoHat3, Rsite1, Rsite2, Rsite3, t1, t2, t3, r2_true)

global mu wEarthVec

mu=3.986e14;

% Velocity and acceleration of site at time t2
RsiteDot2 = cross(wEarthVec, Rsite2);
RsiteDDot2= cross(wEarthVec, RsiteDot2);

% LAPLACE METHOD

% rhoDot(t) at any time
rhoHatDot1 =rhoHatDot_t( t1, t1, t2, t3, rhoHat1, rhoHat2, rhoHat3);
rhoHatDot2 =rhoHatDot_t( t2, t1, t2, t3, rhoHat1, rhoHat2, rhoHat3);
rhoHatDot3 =rhoHatDot_t( t3, t1, t2, t3, rhoHat1, rhoHat2, rhoHat3); 

% rhoDDot(t) at any time
rhoHatDDot1 =rhoHatDot_t( t1, t1, t2, t3, rhoHatDot1, rhoHatDot2, rhoHatDot3);
rhoHatDDot2 =rhoHatDot_t( t2, t1, t2, t3, rhoHatDot1, rhoHatDot2, rhoHatDot3);
rhoHatDDot2 =rhoHatDot_t( t3, t1, t2, t3, rhoHatDot1, rhoHatDot2, rhoHatDot3);

% Determinants of the system of equaitons in slant range
matrixD= [rhoHat2 ; rhoHatDot2 ; rhoHatDDot2]';  detD = 2*det(matrixD);
matrixD1= [rhoHat2 ; rhoHatDot2 ; RsiteDDot2 ]';  detD1 = det(matrixD1);
matrixD2= [rhoHat2 ; rhoHatDot2 ; Rsite2     ]';  detD2 = det(matrixD2);

CC=dot(rhoHat2, Rsite2);

% Coefficients of the 8th order polynomial
a6 = -4*detD1^2/detD^2 + 4*detD1*CC/detD - (norm(Rsite2))^2;
a3 = - 8*mu*detD1*detD2/detD^2 + 4*mu*detD2*CC/detD;
a0 = -4*mu^2*detD2^2/detD^2;

% Roots of the polynomial
Roots_r2 = roots([1 0 a6 0 0 a3 0 0 a0]);

% Positive roots of the polynomial
r2 = posroot(Roots_r2, r2_true);

rho2 = -2*detD1/detD - 2*mu*detD2/(r2^3*detD);

% Finding rhodot
matrixD3= [rhoHat2 ; RsiteDDot2 ; rhoHatDDot2 ]'; detD3 = det(matrixD3);
matrixD4= [rhoHat2 ; Rsite2     ; rhoHatDDot2 ]'; detD4 = det(matrixD4);

rho2dot = - detD3/detD - mu*detD4/(r2^3*detD);

Robject2 = rho2*rhoHat2 + Rsite2;
Vobject2 = rho2dot*rhoHat2 + rho2*rhoHatDot2 + RsiteDot2;
return
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