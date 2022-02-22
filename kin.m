%% A version of the kinematics for the 2RPS - 1 PRS parallel actuated platform

syms theta11 theta22 theta33 p1 p2 p3 h R H L1 L2 L3 alpha phi psi pz;

%%
eq1 = h*sin(phi)*sin(psi) + h*cos(phi) - R + L1*cos(theta11) + p1*cos(theta11) == 0;
eq2 = h*sin(phi)*sin(psi) - h*cos(phi) + R - L3*cos(theta33) - p3*cos(theta33) == 0;
eq3 = h*cos(psi) - p2*cos(alpha) - H*sin(alpha) + L2*cos(theta22 + alpha) - R == 0;
eq4 = pz + h*sin(phi) - L3*sin(theta33) + p3*sin(theta33) + H == 0;
eq5 = pz - h*sin(phi) - L1*sin(theta11) + p1*sin(theta11) + H == 0;
eq6 = pz - h*cos(phi)*sin(psi) + p2*sin(alpha) - L2*sin(theta22 + alpha) - H*cos(alpha) == 0;

%%
assume(p1,'real');
assume(p2,'real');
assume(p3, 'real');
assume(theta11,'real');
assume(theta22, 'real');
assume(theta33, 'real');
sol = solve([eq1,eq2, eq3, eq4, eq5, eq6], [p1, p2, p3, theta11, theta22, theta33]);


%%

%% Forward Kinematics

%Given input p1, p2, p3. Solve for Pz, theta1, theta2.

h = 0.12284;
R = 0.1361;
L1 = 0.10924;
L3 = 0.10924;
H = 0.01908;
L2 = 0.25;
p3 = out.p3.Data;
theta33 = out.theta33.Data;
theta11 = out.theta11.Data;
theta11 = theta11 - 1.56;
p1 = out.p1.Data;
p2 = 100 - out.p2.Data;
alpha = (180-110.74)*3.14/180;
Theta1 = [];
pz = out.Pz.Data;
phi = out.thetay.Data;
psi = out.thetax.Data;
Theta1dup = [];
%h.*cos(psi(i)) - p2.*cos(alpha) - H.*sin(alpha) + L2.*cos(theta22 + alpha)
theta22 = 3.14 - out.theta22.Data;
Theta2 = [];
Pz = [];
for k = 1:length(out.tout)
    Z = p2(k).*cos(alpha) + H.*sin(alpha) - L2.*cos(theta22(k) + alpha) + R;
    %(Z+h/h )t2.^2 + (Z/h - 1) = 0
    p = [((Z/h) + 1) 0 ((Z/h) -1)];
    r = roots(p);
    theta2 = real(2*atan(r));
    Theta2 = [Theta2 theta2(1)];
end
%%
Theta1 = [];
Theta1dup = [];
for j = 1:length(out.tout)
    Z1 = 2*h.*sin(Theta2(j)) + h + R - L3.*cos(theta33(j)) - p3(j).*cos(theta33(j));
    Z2 = h - R + L3*cos(theta33(j)) + p3(j).*cos(theta33(j));
    r1 = roots([Z1 0 -Z2]);
    theta1 = (2*atan(r1));
    Theta1 = [Theta1 theta1(1)];
end
for j = 1:length(out.tout)
    Z1 = 2*h.*sin(Theta2(j)) + h + R - L3.*cos(theta33(j)) - p3(j).*cos(theta33(j));
    Z2 = h - R + L3*cos(theta33(j)) + p3(j).*cos(theta33(j));
    r2 = roots([Z1 0 -Z2]);
    theta1dup = (2*atan(r2));
    Theta1dup = [Theta1dup theta1dup(1)];
end

%%

for g = 1:length(out.tout)
    A1 = L3.*sin(theta33) + p3.*sin(theta33) + H;
    A2 = L1.*sin(theta11) + p1.*sin(theta11) + H;
    A3 = -p2.*sin(alpha) + L2.*sin(theta22 + alpha) + H*cos(alpha);
    A = A1.*A2.*A3;
    Z33 = (A - h.*(sin(Theta1).^2) .*(h.*cos(Theta1).*sin(Theta2)));
    Z31 = -h.*cos(Theta1).*sin(Theta2);
    Z32 = -(h.*sin(Theta1)).^2;
    Z30 = 1;
    r3 = roots([Z30 Z31 Z32 -Z33]);
    pz = (real(r3));
    Pz = [Pz pz];
    Z30*Pz^3 + Z31*Pz^2 + Z32*Pz - Z33 = 0
end

%% Inverse Kinematics
%Given input Pz, theta1 theta2. Solve for p1, p2, p3, theta11, theta22,
%theta33.

Angle11 = [];
for q = 1:length(out.tout)
    Z4 = 2*h.*sin(Theta1dup(q)).*sin(Theta2(q)) + 2*h.*cos(Theta1dup(q)) - 2*R - Pz(q) + h.*sin(Theta1dup(q)) +H;
    Z5 = Pz(q) - h.*sin(Theta1dup(q)) - H;
    T11 = roots([Z4 0 Z5]);
    angle11 = real(2*atan(T11));
    Angle11 = [Angle11 angle11(1)];
end

%%
Angle33 = [];
for w = 1:length(out.tout)
    Z6 = 2*h.*sin(Theta1dup(w)).*sin(Theta2(w)) - 2*h.*cos(Theta1dup(w)) + 2*R + Pz(w) + h.*sin(Theta1dup(w)) - H;
    Z7 = - (Pz(w) + h.*sin(Theta1dup(w)) - H);
    T33 = roots([Z6 0 Z7]);
    angle33 = real(2*atan(T33));
    Angle33 = [Angle33 angle33(1)];
end
%%
Pris1 = [];
for e = 1:length(out.tout)
    Z8 = (h.*sin(Theta1dup(e)).*sin(Theta2(e)) + h.*cos(Theta1dup(e)) - R + L1.*cos(Angle11(e))).^2 + (Pz(e) - h.*sin(Theta1dup(e)) - L1.*sin(Angle11(e)) - H).^2;
    pris1 = roots([1 0 -Z8]);
    Pris1 = [Pris1 pris1(1)];
end

Pris3 = [];
for r = 1:length(out.tout)
    Z9 = (h.*sin(Theta1dup(r)).*sin(Theta2(r)) - h.*cos(Theta1dup(r)) + R - L3.*cos(Angle33(r))).^2 + (Pz(r) + h.*sin(Angle11(r)) - L3.*sin(Angle33(r)) -H).^2 ;
    pris3 = roots([1 0 -Z9]);
    Pris3 = [Pris3 pris3(1)];
end
%%
Angle22 = [];
for t = 1:length(out.tout)
    Z10 = H -Pz(t).*cos(alpha) + h.*cos(Theta1dup(t)).*sin(Theta2(t)).*cos(alpha) - h.*cos(Theta2(t)).*sin(alpha) + R.*sin(alpha);
    Z11 = Z10/(-Z10 - 2*L2);
    T22 = roots([1 0 -Z11]);
    angle22 = (2*atan(T22));
    Angle22 = [Angle22 angle22(1)];
end

%%
Pris2 = [];
for y = 1:length(out.tout)
    Z12 = (h.*cos(Theta2(y)) - H.*sin(alpha) + L2.*cos(Angle22(y) + alpha) - R).^2 + (Pz(y) - h.*cos(Theta1dup(y)).*sin(Theta2(y)) - H.*cos(alpha) - L2.*sin(Angle22(y) + alpha)).^2;
    pris2 = roots([1 0 -Z12]);
    Pris2 = [Pris2 pris2(2)];
end
