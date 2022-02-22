%% Forming all the vectors and kinematic equations for the 2PRU - 1 PRS parallel platform


R = 0.35;
syms beta alpha;
theta = beta;
psi = alpha;
RotMF = [cos(beta) sin(beta)*sin(alpha) sin(beta)*cos(alpha); 0 cos(alpha) -sin(alpha); -sin(beta) cos(beta)*sin(alpha) cos(beta)*cos(alpha)];
a1 = [R; 0; 0];
a2 = [0; R; 0];
a3 = [-R; 0; 0];
A1P = RotMF*a1;
A2P = RotMF*a2;
A3P = RotMF*a3;
 

%%
%kinematics
 A1Px = A1P(1);
 A1Pz = A1P(3);
 A3Px = A3P(1);
 A3Pz = A3P(3);
 A2Py = A2P(2);
 A2Pz = A2P(3);
%%
%%limb1
H = 0.4;
R = 0.35;
syms px py pz H rho1 L theta11 theta12 theta13 A1Px A1Py A1Pz;
P1 = [px; py; pz];
H1matrix = [H; 0; 0];
Rho1matrix = [0; 0; rho1];
L1matrix = [L*cos(theta11); 0; L*sin(theta11)];
Rot1 = [cos(theta12) sin(theta12)*sin(theta13) sin(theta12)*cos(theta13); 0 cos(theta13) -sin(theta13); -sin(theta12) cos(theta12)*sin(theta13) cos(theta12)*cos(theta13)];
AP1mtx = [A1Px; A1Py; A1Pz];
P1 = H1matrix + Rho1matrix + L1matrix + Rot1*AP1mtx;

%%
%limb3
syms px py pz H rho3 L theta31 theta32 theta33 A3Px A3Py A3Pz;
P3 = [px; py; pz];
H3matrix = [-H; 0; 0];
Rho3matrix = [0; 0; rho3];
L3matrix = [L*cos(theta31); 0; L*sin(theta31)];
Rot3 = [cos(theta32) sin(theta32)*sin(theta33) sin(theta32)*cos(theta33); 0 cos(theta33) -sin(theta33); -sin(theta32) cos(theta32)*sin(theta33) cos(theta32)*cos(theta33)];
AP3mtx = [A3Px; A3Py; A3Pz];
P3 = H3matrix + Rho3matrix + L3matrix + Rot3*AP3mtx;

%%
%limb2
syms px py pz H rho2 L theta21 theta22 theta23 theta24 A2Px A2Py A2Pz;
P2 = [px; py; pz];
H2matrix = [H; 0; 0];
Rho2matrix = [0; 0; rho2];
L2matrix = [L*cos(theta21); 0; L*sin(theta21)];
Rot2 = [cos(theta23)*cos(theta24),     cos(theta23)*sin(theta24),   sin(theta23); - cos(theta22)*sin(theta23) - cos(theta24)*sin(theta22)*sin(theta23),   cos(theta22)*cos(theta24) - sin(theta22)*sin(theta23)*sin(theta24), cos(theta23)*sin(theta22) ; sin(theta22)*sin(theta24) - cos(theta22)*cos(theta24)*sin(theta23), - cos(theta24)*sin(theta22) - cos(theta22)*sin(theta23)*sin(theta24), cos(theta22)*cos(theta23)];
AP2mtx = [A2Px; A2Py; A2Pz];
P2 = H2matrix + Rho2matrix + L2matrix + Rot2*AP2mtx;

 %% Forward Kinematics
 %value for beta
 beta = atan2(-A2Py, A2Pz);
 calpha1 = (L*cos(theta11) - H - A1Px*cos(beta))/A1Pz*sin(beta);
 calpha2 = (H - L*cos(theta31) - A3Px*cos(beta))/(A3Pz*sin(beta));
 
 salpha = (A2Py*calpha1 + H - L*cos(theta21))/(A2Pz*cos(beta) + A2Py*sin(beta));
 
 alpha = atan2(salpha, calpha1);
 %%
 Pz1 = rho1 + L*sin(theta21) + A2Pz*cos(alpha)*cos(beta) - A2Py*sin(alpha) - A2Py*cos(alpha)*sin(beta);
 Pz2 = rho2 + L*sin(theta11) - A1Px *sin(beta) + A1Pz*cos(beta)*cos(alpha);
 Pz3 = rho3 + L*sin(theta31) -A3Px*sin(beta) + A3Pz*cos(beta)*cos(alpha);


%% Inverse
 %theta = beta, psi = alpha
 R = 0.25;
 H = 0.40;
 z = 0.02;
 L = 0.4;
 Rho1 = [];
 Rho2 = [];
 Rho3 = [];
 theta = [-0.8722:0.01:0.8722];
 psi = [-0.8722:0.01:0.8722];
 for k = 1:length(theta)
     rho1 = z - R*sin(theta(k)) - sqrt(L.^2 - (R*(cos(theta(k)) - sin(theta(k))*sin(psi(k))) -H).^2);
     rho2 = z + R*cos(theta(k))*sin(psi(k)) - sqrt(L.^2 - (R*cos(psi(k)) -H).^2);
     rho3 = z + R*sin(theta(k))- sqrt(L.^2 - (-R*(cos(theta(k)) + sin(theta(k))*sin(psi(k))) + H).^2);
     Rho1 = [Rho1 , rho1];
     Rho2 = [Rho2 ,rho2];
     Rho3 = [Rho3, rho3];
 end