%% RPS Ankle rehabilitation Kinematics (Forward Kinematics using a Quaternion Approach)


syms theta psi phi X Y Z eps2 eps3 eps1 a1 a2 a3 b1 b2 b3;
%Constraints
syms x0 x1 x2 x3;
eps2 = 1.046;
eps3 = -1.046;
a3 = 0.17785;
b3 = 0.07684;

Rot = [x0.^2 + x1.^2 + x2.^2 + x3.^2, 0,0,0; X, x0.^2 + x1.^2 - x2.^2 - x3.^2, 2*x1*x2 - 2*x0*x3, 2*x1*x3 + 2*x0*x2; Y, 2*x1*x2 + 2*x0*x3, x0.^2 - x1.^2 + x2.^2 - x3.^2, 2*x3*x2 - 2*x0*x1; Z, 2*x3*x1 - 2*x0*x2, 2*x3*x2 + 2*x0*x1, x0.^2 - x1.^2 - x2.^2 + x3.^2];   

%% Design 2C
a1 = -a3*cos(eps2 - eps3)/cos(eps2);
a2 = a3*cos(eps3)/cos(eps2);
b1 = -b3*cos(eps2 - eps3)/cos(eps2);
b2 = b3*cos(eps3)/cos(eps2);

%% Forming constraint equations - 2
%%Rotation matrix when x0 = x(1), x1 = x(2), x2 = x(3), x3 =x(4), X =
%%x(5),Y = x(6), Z = x(7)
Rotx10 = [x0.^2 + x1.^2 + x2.^2 + x3.^2, 0,0,0; X, x0.^2 + x1.^2 - x2.^2 - x3.^2, 2*x1*x2 - 2*x0*x3, 2*x1*x3 + 2*x0*x2; Y, 2*x1*x2 + 2*x0*x3, x0.^2 - x1.^2 + x2.^2 - x3.^2, 2*x3*x2 - 2*x0*x1; Z, 2*x3*x1 - 2*x0*x2, 2*x3*x2 + 2*x0*x1, x0.^2 - x1.^2 - x2.^2 + x3.^2];   

 pvalue = 0.1201;
 ri = 0.06174;
 ri = ri +pvalue; 
 r0A1 = [1; a1; 0; 0];
 r0A2 = [1; cos(eps2)*a2; sin(eps2)*a2; 0];
 r0A3 = [1; cos(eps3)*a3; sin(eps3)*a3; 0];
 
 r1B1 = [1; b1; 0; 0];
 r1B2 = [1; cos(eps2)*b2; sin(eps2)*b2; 0];
 r1B3 = [1; cos(eps3)*b3; sin(eps3)*b3; 0];
 
 r0B1 = Rotx10*r1B1;
 r0B2 = Rotx10*r1B2;
 r0B3 = Rotx10*r1B3;
 eq1 = cos(eps2)*Y - 2*cos(eps2 - eps3)*(x0*x3 + x1*x2)*b3 == 0;
 eq2 = -sin(eps2)*cos(eps2)*X + cos(eps2).^2 * Y + (- 2*cos(eps3)*cos(eps2)*sin(eps2)*x1.^2 -sin(eps2)*cos(eps2)*x2.^2 - 2*cos(eps2).^2*x2*x1 - x0*x3 + x1*x2)*b3 == 0 ;
 eq3 = -sin(eps3)*X + cos(eps3)* Y + (- 2*cos(eps3)*sin(eps3)*x1.^2 +2*sin(eps3)*cos(eps3)*x2.^2 + 4*cos(eps3).^2*x2*x1 + 2*x0*x3 -2*x1*x2)*b3 == 0 ;

 eq4 = (r0B1(1) - r0A1(1)).^2 + (r0B1(2) - r0A1(2)).^2 + (r0B1(3) - r0A1(3)).^2 + (r0B1(4) - r0A1(4)).^2 - ri.^2 == 0; 
 
eq5 = (r0B2(1) - r0A2(1)).^2 + (r0B2(2) - r0A2(2)).^2 + (r0B2(3) - r0A2(3)).^2 + (r0B2(4) - r0A2(4)).^2 - ri.^2 == 0;
 
eq6 = (r0B3(1) - r0A3(1)).^2 + (r0B3(2) - r0A3(2)).^2 + (r0B3(3) - r0A1(3)).^2 + (r0B3(4) - r0A3(4)).^2 - ri.^2 == 0;

%% Non linear equation solver.... x0 = x(1), x2 = x(2), x3 = x(3), X = x(4), Y = x(5), Z = x(6)
fun = @root2d;
x0 = [0.4, 0.4, 0.4,0.4,  0.12, 0.12, 0.12];
options = optimoptions('fsolve','Display','iter');
options.StepTolerance = 1.0e-6;
options.FunctionTolerance = 1.00;
sol  = fsolve(fun, x0, options);

%% Rotation matrix when x2 = 0
x2 = 0;
Rotx20 = [x0.^2 + x1.^2 + x2.^2 + x3.^2, 0,0,0; X, x0.^2 + x1.^2 - x2.^2 - x3.^2, 2*x1*x2 - 2*x0*x3, 2*x1*x3 + 2*x0*x2; Y, 2*x1*x2 + 2*x0*x3, x0.^2 - x1.^2 + x2.^2 - x3.^2, 2*x3*x2 - 2*x0*x1; Z, 2*x3*x1 - 2*x0*x2, 2*x3*x2 + 2*x0*x1, x0.^2 - x1.^2 - x2.^2 + x3.^2];   

 pvalue = 0.1201;
 ri = 0.06174;
 ri = ri +pvalue; 
 r0A1 = [1; a1; 0; 0];
 r0A2 = [1; cos(eps2)*a2; sin(eps2)*a2; 0];
 r0A3 = [1; cos(eps3)*a3; sin(eps3)*a3; 0];
 
 r1B1 = [1; b1; 0; 0];
 r1B2 = [1; cos(eps2)*b2; sin(eps2)*b2; 0];
 r1B3 = [1; cos(eps3)*b3; sin(eps3)*b3; 0];
 
 r0B1 = Rotx20*r1B1;
 r0B2 = Rotx20*r1B2;
 r0B3 = Rotx20*r1B3;
 
 

eq4 = (r0B1(1) - r0A1(1)).^2 + (r0B1(2) - r0A1(2)).^2 + (r0B1(3) - r0A1(3)).^2 + (r0B1(4) - r0A1(4)).^2 - ri.^2 == 0; 
 
eq5 = (r0B2(1) - r0A2(1)).^2 + (r0B2(2) - r0A2(2)).^2 + (r0B2(3) - r0A2(3)).^2 + (r0B2(4) - r0A2(4)).^2 - ri.^2 == 0;
 
eq6 = (r0B3(1) - r0A3(1)).^2 + (r0B3(2) - r0A3(2)).^2 + (r0B3(3) - r0A1(3)).^2 + (r0B3(4) - r0A3(4)).^2 - ri.^2 == 0;
%% 
sol = solve(eq1,eq2,eq3,eq4,eq5,eq6,[x0, x1, x2, x3,X,Y]);

%% Solution for theta and phi when x2 = 0

x0 = cos(theta/2)*cos(phi);
x1 = sin(theta/2);
x3 = cos(theta/2)*sin(phi);

%%  Solution for theta and phi when x1 = 0
x0 = cos(theta/2) *sin(phi);
x2 = sin(theta/2);
x3 = -cos(theta/2)*cos(phi);
%%

%%
