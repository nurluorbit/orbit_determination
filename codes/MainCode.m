% ------------------------------------------------------------
% Gauss's and Laplace's Methods for
% Preliminary Orbit Determination
% pi - 3.1415926...
% mu - gravitational parameter (mˆ3/sˆ2)
% coe - orbital elements [h e RA incl w TA a]
% where h = angular momentum (mˆ2/s)
% e = eccentricity
% RA = right ascension of the ascending node (rad)
% incl = orbit inclination (rad)
% w = argument of perigee (rad)
% TA = true anomaly (rad)
% a = semimajor axis (m)
% r - position vector (m) in geocentric equatorial frame
% v - velocity vector (m/s) in geocentric equatorial frame
% f - The oblateness, or flattening (unitless)
% H - Sea level /elevation/
% Rsoi - The radius of the sphere of influence
% TAsoi - Maximum true anomaly (rad)
% ------------------------------------------------------------
close all; clc; clear;
% Coefficients
global mu
global wEarthVec 
format long

% Coefficients=[];
deg = pi/180;
mu = 3.98597913*10^14; % [m^3/s^2]
Re = 6.3781366*10^6; % [m]
% Earth's rotational speed and its vector
wEarth = 2*pi/86164.1;
wEarthVec = wEarth*[0 0 1];
f = 0;
H = 0; % [m]

% Given orbital elements' values
h = 8*10^10; % [m^2/s]
h = 7.7783969e10;
ecc = 1.2;
RA = 60; % [deg]
incl = 90; %80 % [deg]
w = 40; % [deg]
if ecc < 1
    asemi = h^2 / mu / (1 - ecc^2);
    Tper = 2*pi * sqrt(asemi^3 / mu);
else
    asemi = 0;
    Tper = 0;
end

%%Convert degrees to radians
RArad = RA*deg; % [rad]
inclrad = incl*deg; % [rad]
wrad = w*deg; % [rad] 

% Longitude and latitude
Lon = 30; % [deg]
Lat = 89; % [deg]

% Date and time
year = 2018;
month = 3;
day = 16;
hour = 9;
minute = 32;
second = 27;

% Universal time
ut = hour + minute/60 + second/3600;

% Equation 5.48 from Curtis
j0 = J0(year, month, day);

% Equation 5.47 from Curtis
jd = j0 + ut/24;

% Curtis equation 5.54
Rphi = Re / sqrt(1 - (2*f - f^2)*sin(Lat)^2);

% Curtis equation 5.55b
Rc = Rphi + H;

% Curtis equation 5.55b
Rs = (1-f)^2*Rphi + H;

Error = 0;

taDiff=500; taLast=1e6; taFirst=-taLast;
if ecc < 1
    taFirst = -taDiff;
    taLast = Tper+taDiff;
end
Gozlem=[];
for t=taFirst:taDiff:taLast
    
% True Anomaly
TA = TimeSec2TArad(h, ecc, t, mu, 1e-7);
    
% Julian days
jdi = jd + t/(60*60*24);
    
% Universal time
 uti = (jdi - j0)*24;
    
% Local sidereal time
lsti = LST(j0, uti, Lon);
    
% The observer’s position vector [m]
Rsitei = [Rc*cosd(Lat)*cosd(lsti)  Rc*cosd(Lat)*sind(lsti) Rs*sind(Lat)];
    
% Orbital elements
coe_i = [h, ecc, RArad, inclrad, wrad, TA];
    
% Calculation the position vectors from orbital elements [m]-[m/s]
[ri, vi] = sv_from_coe(coe_i);
    
% The magnitude of the given position vector
ri_norm = norm(ri);
    
% The slant range / The topocentric direction cosine vector
rtopoi = ri - Rsitei;
    
% calculation topocentric right ascension and declination for Gauss method
[alfatopoi, deltatopoi] = RADEC_from_r(rtopoi);
    
alfatopoi = (1 + (rand - 0.5)*2*Error)*alfatopoi;
deltatopoi = (1 + (rand - 0.5)*2*Error)*deltatopoi;
    
% The slant range (rtopoi)
rho = [cosd(deltatopoi)*cosd(alfatopoi), cosd(deltatopoi)*sind(alfatopoi), sind(deltatopoi)];
    
% Observation's variables
Gozlem=[Gozlem; t, TA, lsti, alfatopoi, deltatopoi, Rsitei, rho, ri_norm];

TA    
end

sizeGozlem=size(Gozlem);

for i=1:sizeGozlem(1)-2
    [sizeGozlem i]
    
% Separate observation's variables into three parts
Rho1 = Gozlem(i, 9:11);
Rho2 = Gozlem(i+1, 9:11);
Rho3 = Gozlem(i+2, 9:11);
    
R1(i,:) = Gozlem(i, 6:8);
R2(i,:) = Gozlem(i+1, 6:8);
R3(i,:) = Gozlem(i+2, 6:8);
    
t1(i) = Gozlem(i, 1);
t2(i) = Gozlem(i+1, 1);
t3(i) = Gozlem(i+2, 1);
    
% Find the state vectors without iterative improvement using Gauss
[rG(i,:), vG(i,:)] = Gauss(Rho1, Rho2, Rho3, R1(i,:), R2(i,:),R3(i,:), t1(i), t2(i), t3(i), Gozlem(i+1, 12));
    
coe_calcG(i, 1) = t2(i);
% Calculate orbital elements from the state vectors without iterative improvement
coe_calcG(i, 2:8) = coe_from_sv(rG(i,:), vG(i,:));
end

for i=1:sizeGozlem(1)-2
    [sizeGozlem i]
    
% Separate observation's variables into three parts
Rho1 = Gozlem(i, 9:11);
Rho2 = Gozlem(i+1, 9:11);
Rho3 = Gozlem(i+2, 9:11);
    
R1(i,:) = Gozlem(i, 6:8);
R2(i,:) = Gozlem(i+1, 6:8);
R3(i,:) = Gozlem(i+2, 6:8);
    
t1(i) = Gozlem(i, 1);
t2(i) = Gozlem(i+1, 1);
t3(i) = Gozlem(i+2, 1);
    
% Find the state vectors without iterative improvement using Laplace
[rL(i,:), vL(i,:)] = Laplace(Rho1, Rho2, Rho3, R1(i,:), R2(i,:), R3(i,:), t1(i), t2(i), t3(i), Gozlem(i+1, 12));

coe_calcL(i, 1) = t2(i);
% Calculate orbital elements from the state vectors without iterative improvement
coe_calcL(i, 2:8) = coe_from_sv(rL(i,:), vL(i,:));
end

%% Rsite and orbital plane check
for i=1:sizeGozlem(1)-2
    kisiiL(i) = pi/2 - acos(dot(cross(rL(i,:)./norm(rL(i,:)),vL(i,:)./norm(vL(i,:))), R2(i,:)./norm(R2(i,:))));
end
kisiiL = kisiiL';

%% Rsite and orbital plane check
for i=1:sizeGozlem(1)-2
    kisiiG(i) = pi/2 - acos(dot(cross(rG(i,:)./norm(rG(i,:)),vG(i,:)./norm(vG(i,:))), R2(i,:)./norm(R2(i,:))));
end
kisiiG = kisiiG';
%% Figure 1
if ecc<1
EndAngle=2*pi;
else
EndAngle=acos(-1/ecc)-1e-3;
end;
TAg = -EndAngle:0.01:EndAngle;
r = h^2 ./ mu ./ (1 + ecc .* cos(TAg));
x = r .* cos(TAg);
y = r .* sin(TAg);
figure(1)
plot(x , y, 'k');
r = coe_calcG(:,2).^2 / mu ./ (1 + coe_calcG(:,3).*cos(coe_calcG(:,7)));
xG = r .* cos(coe_calcG(:,7));
yG = r .* sin(coe_calcG(:,7));
hold on;
plot(xG, yG, 'b');
r = coe_calcL(:,2).^2 / mu ./ (1 + coe_calcL(:,3).*cos(coe_calcL(:,7)));
xL = r .* cos(coe_calcL(:,7));
yL = r .* sin(coe_calcL(:,7));
hold on;
plot(xL, yL, 'r');
angle=0:0.1:2.1*pi;
hold on;
plot(Re*cos(angle),Re*sin(angle), 'm');
axis square;
legend('Given Trajectory','Gauss','Laplace','The Earth ','Location','best');
if ecc < 1
title('Trajectory - Elliptical Case');
rapo = h^2/mu /(1-ecc);
xlim([-rapo*1.1 rapo*1.1])
ylim([-rapo rapo])
elseif ecc == 1
title('Trajectory - Parabolic Case');
elseif ecc > 1
title('Trajectory - Hyperbolic Case');
end
grid on;
grid minor;
%% Figure 2
figure(2)
semilogy(coe_calcG(:,1),abs(h-coe_calcG(:,2))/h, 'b');
hold on;
semilogy(coe_calcL(:,1),abs(h-coe_calcL(:,2))/h, 'r');
hold on;
legend('Gauss','Laplace','Location','best');
if ecc < 1
title('Relative Error in Angular Momentum - Elliptic Case');
elseif ecc == 1
title('Relative Error in Angular Momentum - Parabolic Case');
elseif ecc > 1
title('Relative Error in Angular Momentum - Hyperbolic Case');
end
xlabel('Time [sec]','FontSize',12,'FontWeight','bold','Color','k');
%ylabel('h ','FontSize',12,'FontWeight','bold','Color','k');
grid on;
grid minor;
axis square;
%% Figure 3
figure(3)
hold on;
semilogy(coe_calcG(:,1),abs(ecc-coe_calcG(:,3))/ecc, 'b');
semilogy(coe_calcL(:,1),abs(ecc-coe_calcL(:,3))/ecc, 'r');
legend('Gauss','Laplace','Location','best');
if ecc < 1
title('Relative Error in Eccentricity - Elliptic Case');
elseif ecc == 1
title('Relative Error in Eccentricity - Parabolic Case');
elseif ecc > 1
title('Relative Error in Eccentricity - Hyperbolic Case');
end
xlabel('Time [sec]','FontSize',12,'FontWeight','bold','Color','k');
%ylabel('e','FontSize',12,'FontWeight','bold','Color','k');
grid on;
grid minor;
axis square;
%% Figure 4
figure(4)
semilogy(coe_calcG(:,1),abs(RArad-coe_calcG(:,4))/RArad, 'b');
hold on;
semilogy(coe_calcL(:,1),abs(RArad-coe_calcL(:,4))/RArad, 'r');
hold on;
legend('Gauss','Laplace','Location','best');
if ecc < 1
title('Relative Error in Right Ascension of the Ascending Node - Elliptic Case');
elseif ecc == 1
title('Relative Error in Right Ascension of the Ascending Node - Parabolic Case');
elseif ecc > 1
title('Relative Error in Right Ascension of the Ascending Node - Hyperbolic Case');
end
xlabel('Time [sec]','FontSize',12,'FontWeight','bold','Color','k');
%ylabel('RA [rad]','FontSize',12,'FontWeight','bold','Color','k');
grid on;
grid minor;
axis square;
%% Figure 5 
figure(5)
semilogy(coe_calcG(:,1),abs(inclrad-coe_calcG(:,5))/inclrad, 'b');
hold on;
semilogy(coe_calcL(:,1),abs(inclrad-coe_calcL(:,5))/inclrad, 'r');
hold on;
legend('Gauss','Laplace','Location','best');
if ecc < 1
title('Relative Error in Orbit Inclination - Elliptic Case');
elseif ecc == 1
title('Relative Error in Orbit Inclination - Parabolic Case');
elseif ecc > 1
title('Relative Error in Orbit Inclination - Hyperbolic Case');
end
xlabel('Time [sec]','FontSize',12,'FontWeight','bold','Color','k');
%ylabel('incl [rad]','FontSize',12,'FontWeight','bold','Color','k');
grid on;
grid minor;
axis square;
%% Figure 6  
figure(6)
semilogy(coe_calcG(:,1),abs(wrad-coe_calcG(:,6))/wrad, 'b');
hold on;
semilogy(coe_calcL(:,1),abs(wrad-coe_calcL(:,6))/wrad, 'r');
hold on;
legend('Gauss','Laplace', 'Given Orbit', 'Location','best');
if ecc < 1
title('Relative Error in Argument of Perigee - Elliptic Case');
elseif ecc == 1
title('Relative Error in Argument of Perigee - Parabolic Case');
elseif ecc > 1
title('Relative Error in Argument of Perigee - Hyperbolic Case');
end
xlabel('Time [sec]','FontSize',12,'FontWeight','bold','Color','k');
% ylabel('w [rad]','FontSize',12,'FontWeight','bold','Color','k');
grid on;
grid minor;
axis square;
%  %% Figure 7
% figure(7)
% hold on;
% plot(coe_calcL(:,1),coe_calcL(:,3), 'r');
% plot(coe_calcG(:,1),coe_calcG(:,3), 'b');
% plot(coe_calcG(:,1), kisiiG, 'm');
% plot(coe_calcL(:,1), kisiiL, 'k');
% legend('Laplace', 'Gauss', 'Location','best');
% title('Failing Moments of the Methods ');
% xlabel('Time [sec]','FontSize',12,'FontWeight','bold','Color','k');
% ylabel('Eccentricity, e ','FontSize',12,'FontWeight','bold','Color','k');
% grid on;
% grid minor;
% axis square;