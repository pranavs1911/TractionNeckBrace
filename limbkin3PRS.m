%% Numerical Solutions for the 3PRS mechanism
%Kinematics

syms theta phi psi;
syms r1x r1y r1z r2x r2y r2z r3x r3y r3z;
syms x y z;
syms rp alpha beta;
smalla1dash = [rp; 0;0];
smalla2dash = [rp*cos(alpha); rp*sin(alpha); 0];
smalla3dash = [rp*cos(beta); rp*sin(beta); 0];
T = [cos(theta)*cos(phi) + sin(psi)*sin(theta)*sin(phi), -cos(theta)*sin(phi) + sin(psi)*sin(theta)*cos(phi) , cos(psi)*sin(theta); cos(psi)*sin(phi), cos(psi)*cos(phi), -sin(psi); -sin(theta)*cos(phi) + sin(psi)*cos(theta)*sin(phi), sin(theta)*sin(phi) + sin(psi)*cos(theta)*cos(phi) , cos(psi)*cos(theta)];
P = [x; y; z];
A1 = T*smalla1dash;
A2 = T*smalla2dash;
A3 = T*smalla3dash;
r1 = P + A1;
r2 = P + A2;
r3 = P + A3;
r1x = r1(1);
r1y = r1(2);
r1z = r1(3);
r2x = r2(1);
r2y = r2(2);
r2z = r2(3);

r3x = r3(1);
r3y = r3(2);
r3z = r3(3);

%% Constraints
r2y = r2x*tan(alpha);
r3y = r3x*tan(beta);

%% Calculation of phi
%psi - X axis, theta - Yaxis, 
alpha = 90*3.14/180;
beta = 270*3.14/180;

theta = out.thetax.Data;
psi = out.thetaz.Data;
phil = [];
for i = 1:length(psi)
    A = cos(alpha) - cos(beta);
    B = sin(alpha) - sin(beta);
    C = ((cos(beta)-1)/tan(beta)) - ((cos(alpha) - 1)/tan(alpha));
    R = A*(cos(theta(i)) - cos(psi(i))) + B*(cos(theta(i))).*sin(psi(i));
    S = A*sin(theta(i)).*sin(psi(i)) - B*cos(theta(i)) + C*cos(psi(i));
    phi = atan2(-R,S);
    phil = [phil phi];
end




%% Calculation of x and y
rp = 0.055;
theta = out.thetax.Data;
psi = out.thetaz.Data;
Xmtx = [];
Ymtx =[];
for l = 1:length(theta)
    xi = -rp*(cos(theta(l))*cos(phil(l)) + sin(theta(l))*sin(psi(l))*sin(phil(l)))*cos(alpha) + (-rp)*(-cos(theta(l))*sin(phil(l)) + sin(psi(l))*sin(theta(l))*cos(phil(l)))*sin(alpha);  
    partx = ((rp/tan(alpha)) *(cos(psi(l))*sin(phil(l)) *(cos(alpha) - 1)  + cos(psi(l))*cos(phil(l))*sin(alpha)));
    y = -cos(psi(l))*sin(phil(l))*rp;
    x = xi + partx;
    Xmtx = [Xmtx , x];
    Ymtx = [Ymtx, y];
end

%% Inverse Kinematics analysis Actuator 1
rp = 0.055;
theta = out.thetax.Data; 
psi = out.thetay.Data;
phil = out.phi.Data;
l = 0.17088;
y = Ymtx;
x = Xmtx;
z = 0.02;
Act1 = [];
i = 1;
%Calculation of b1, b2 ,b3
while i <= length(theta)
    r1x = x(i) + rp*(cos(phil(i))*cos(theta(i)) + sin(phil(i))*sin(psi(i))*sin(theta(i)));
    r1y = y(i) + rp*cos(psi(i))*sin(phil(i));
    r1z = z - rp*(cos(phil(i))*sin(theta(i)) - cos(theta(i))*sin(phil(i))*sin(psi(i)));
    A1 = 1;
    B1 = -2*r1x;
    C1 = r1x.^2 + r1y.^2 + r1z.^2 - l.^2 ;
    b1 = roots([A1 B1 C1]);
    b1 = b1 (b1>0);
    Act1 = [Act1, b1];
    i = i + 1;
end

%%  Inverse Kinematics analysis Actuator 2

rp = 0.055;
theta = out.thetax.Data; 
psi = out.thetay.Data;
phil = out.phi.Data;
l = 0.158;
y = Ymtx;
x = Xmtx;
z = 0.02;
Act2 = [];
j = 1;
%Calculation of b1, b2 ,b3
while j<=length(theta)
    r2x = x(j) + rp*cos(alpha)*(cos(phil(j))*cos(theta(j)) + sin(phil(j))*sin(psi(j))*sin(theta(j))) - rp*sin(alpha)*(cos(theta(j))*sin(phil(j)) - cos(phil(j))*sin(psi(j))*sin(theta(j)));
    r2y = y(j) + rp*cos(alpha)*cos(psi(j))*sin(phil(j)) + rp*cos(phil(j))*cos(psi(j))*sin(alpha);
    r2z = z - rp*cos(alpha)*(cos(phil(j))*sin(theta(j)) - cos(theta(j))*sin(phil(j))*sin(psi(j))) + rp*sin(alpha)*(sin(phil(j))*sin(theta(j)) + cos(phil(j))*cos(theta(j))*sin(psi(j)));
    A2 = 1;
    B2 = -2*r2x*cos(alpha) - 2*r2y*sin(alpha);
    C2 = r2x.^2 + r2y.^2 + r2z.^2 - l.^2 ;
    b2 = roots([A2 B2 C2]);
    b2 = b2 (b2>0);
    Act2 = [Act2, b2];
    j = j+1;
end

%% Inverse Kinematics analysis Actuator 3

rp = 0.055;
theta = out.theta.Data;  
psi = out.psi.Data; 
phil = out.phi.Data;
l = 0.158;
y = Ymtx;
x = Xmtx;
z = out.z.Data;
Act3 = [];
m = 1;
%Calculation of b1, b2 ,b3
while m <= length(theta)
    r3x =  x(m) + rp*cos(beta)*(cos(phil(m))*cos(theta(m)) + sin(phil(m))*sin(psi(m))*sin(theta(m))) - rp*sin(beta)*(cos(theta(m))*sin(phil(m)) - cos(phil(m))*sin(psi(m))*sin(theta(m)));
    r3y = y(m) + rp*cos(beta)*cos(psi(m))*sin(phil(m)) + rp*cos(phil(m))*cos(psi(m))*sin(beta);
    r3z = z(m) - rp*cos(beta)*(cos(phil(m))*sin(theta(m)) - cos(theta(m))*sin(phil(m))*sin(psi(m))) + rp*sin(beta)*(sin(phil(m))*sin(theta(m)) + cos(phil(m))*cos(theta(m))*sin(psi(m)));
    A3 = 1;
    B3 = -2*r3x*cos(beta) - 2*r3y*sin(beta);
    C3 = r3x.^2 + r3y.^2 + r3z.^2 - l.^2 ;
    b3 = roots([A3 B3 C3]);
    b3 = b3(1);
    Act3 = [Act3, b3];
    m = m+1;
end

