%% Kinematic Solutions for the 2RPS - 1PRS parallel platform
h1 = 0.061455;
g1 = 0.11668;

%% Analytical Solver
syms h1 g1 ux uy uz vx vy vz wx wy wz theta phi psi Pxi Pyi Pzi;
OB1 = [-g1; 0; 0];
OB3 = [0; -g1; 0];
OB2 = [g1; 0; 0];
PS1 = [-h1; 0; 0];
PS2 = [h1; 0; 0];
PS3 = [0; -h1; 0];

%%
phi = 0;
ux = cos(theta)*cos(phi) + sin(psi)*sin(theta)*sin(phi);
uy = cos(theta)*sin(phi) + sin(psi)*sin(theta)*cos(phi);
uz = sin(theta)*cos(psi);
vx = cos(psi)*sin(phi);
vy = cos(psi)*cos(phi);
vz = -sin(psi);
wx = -sin(theta)*cos(phi) + sin(psi)*cos(theta)*sin(phi);
wy = sin(phi)*sin(theta) + sin(psi)*cos(theta)*cos(phi);
wz = cos(theta)*cos(psi);

%%
ORP = [ux, uy, uz;vx, vy, vz;wx, wy, wz];
OP = [Pxi; Pyi; Pzi];
Qvec1 = OP + ORP*PS1;
Qvec2 = OP + ORP*PS2;
Qvec3 = OP + ORP*PS3;

%% Inverse Kinematics:
R1 = 0.13375;
Rev = 0.02394;
d = 67.02*3.14/180;
d1 = 22.98*3.14/180;  
Link3 = [out.Link3x; out.Link3y; out.Link3z];
L1 = Qvec1 - [-R1;0;0] - [0;0;Rev];
L2 = Qvec2 - [R1;0;0] - [0;0;Rev];
L3 = Qvec3 - [0;-0.035;0] - [0; -Rev*cos(d1); Rev*sin(d1)] - Link3;

%% Forward Kinematics
syms theta11 theta21 theta31 L1 L2 L3
R1 = 0.13375;
Rev = 0.02394;
d = 67.02*3.14/180;
d1 = 22.98*3.14/180; 
Link3 = 0.234;
Qvec1dup = [L1*cos(theta11); 0; L1*sin(theta11)]+ [-R1;0;0] + [0;0;Rev];
Qvec2dup = [-L2*cos(theta21); 0; L2*sin(theta21)] + [R1;0;0] + [0;0;Rev];
Qvec3dup = [0;-L3*cos(67.02);-L3*sin(67.02)] + [0;-0.035;0] + [0; -Rev*cos(d1); Rev*sin(d1)] +[0; Link3*cos(theta31); Link3*sin(theta31)];

%% Equations 
%theta
ctheta = ((2*107/800) - (L2*cos(theta21)) - (L1*cos(theta11)))/(2*h1);
stheta = (L1*sin(theta11) - L2*sin(theta21))/(2*h1);
theta = atan2(stheta,ctheta);

%% psi
cpsi = (0.0570 - 0.5005*L3 - 0.234*cos(theta31))/h1;
spsi = (L1*cos(theta11) - (107/800) + h1*cos(theta))/(h1*sin(theta));
psi = atan2(spsi,cpsi);

%% Pzi
Pzi = (L1*sin(theta11) + (1197/50000) + L2*sin(theta21) + (1197/50000))/2;

