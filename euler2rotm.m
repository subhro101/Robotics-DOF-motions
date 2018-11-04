function[ R,D] = euler2rotm(euler)
% assume euler = [phi, theta, psi]
% see function rotm2euler()
%find derivatives here
phi   = euler(1);
theta = euler(2);
psi   = euler(3);
%%drew rooks 
r11 = cos(theta) * cos(phi);
d11= -cos(phi)*sin(theta)-cos(theta)*sin(phi);
r12 = sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi);
d12 = cos(theta)*cos(phi)*sin(psi)-sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)-cos(phi)*cos(psi)+sin(phi)*sin(psi);
r13 = cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi);
d13= cos(theta)*cos(phi)*cos(psi)-cos(psi)*sin(theta)*sin(phi)-cos(phi)*sin(theta)*sin(psi)+ cos(phi)*sin(psi)+cos(psi)*sin(phi);


r21 = cos(theta) * sin(phi);
d21 = -sin(theta)*sin(phi)+cos(theta)*cos(phi);
r22 = sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi);
d22 = cos(theta)*sin(phi)*sin(psi)+cos(phi)*sin(theta)*sin(psi)+cos(psi)*sin(theta)*sin(phi)-cos(psi)*sin(phi)-cos(phi)*sin(psi);
r23 = cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi);
d23 = cos(theta)*cos(psi)*sin(phi)+cos(theta)*cos(psi)*sin(theta)-sin(theta)*sin(phi)*sin(psi)+sin(phi)*sin(psi)-cos(phi)*cos(psi);

r31 = -sin(theta);
d31= -cos(theta);
r32 = sin(psi) * cos(theta);
d32= -sin(theta)*sin(psi)+cos(theta)*cos(psi);
r33 = cos(psi) * cos(theta);
d33 = -cos(psi)*sin(theta)-cos(theta)*sin(psi);

R = [r11, r12, r13;
     r21, r22, r23;
     r31, r32, r33];
D = [d11, d12, d13;
     d21, d22, d23;
     d31, d32, d33];

end