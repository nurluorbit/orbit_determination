function [alfatopo, deltatopo] = RADEC_from_r(r)

% Curtis Orbital Mechanics page 155 (pdf:172)


len_r = norm(r);
% Direction cosines of r
u_r_x = r(1)/len_r;
u_r_y = r(2)/len_r;
u_r_z = r(3)/len_r;
% Declination
deltatopo = asin(u_r_z)*180/pi;
% Right ascension:
if (u_r_y >0)
    alfatopo = acos(u_r_x/cosd(deltatopo))*180/pi;
else
    alfatopo = 360 - acos(u_r_x/cosd(deltatopo))*180/pi;
end
% fprintf('Right ascension = %4.2f [deg] \n',alfa);
% fprintf('Declination = %4.2f [deg] \n',delta);