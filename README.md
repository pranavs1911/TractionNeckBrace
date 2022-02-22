# TractionNeckBrace
Contains all the MATLAB and the simulation files I developed for simulating the initial prototype of the 2RPS - 1 PRS Traction Neck Brace. 

The initial project involved a research of the numerical solution of the forward/ inverse kinematics of the 2RPS - 1PRS parallel platform. 
Various parallel platforms such as 2PRU - 1PRS, 3PRS, 3RPS parallel platforms were analyzed and their kinematic equations were analyzed. 

The 2RPS - 1PRS platform was decided as the ideal structure for the traction neck brace after comparing the degrees of freedom and the kinematics of 
the other parallel actuated platforms mentioned before. 

## Analytical Solution.m 
Developed closed form solutions (forward/inverse kinematics) for the 2RPS - 1PRS parallel actuated platform. 

## kin.m
Contains a version of the analytical kinematics solution for the 2RPS -1PRS parallel platform where I developed the vector matrices
in each limb and analyzed the kinematics.

## limbkin.m
Contains the vector formulas and the kinematic equations for the 2PRU - 1PRS parallel actuated platform. 

## limbkin3PRS.m
Analytic solution for the 3PRS parallel actuated platform. 

## limbkinRPS.m
Contains a derivation of the kinematics for 3RPS ankle rehabilitation device where the closed form solutions are derived using a quaternion approach. 

Reference: Nurahmi, Latifah & Caro, Stéphane & Solichin, Mochamad. (2019). A novel ankle rehabilitation device based on a reconfigurable 3-RPS parallel manipulator. Mechanism and Machine Theory. 134. 135-150. 10.1016/j.mechmachtheory.2018.12.017. 


## Initial Simulation of the draft design on Simulink. 

The 2RPS - 1PRS platform is designed such that it provides three independent motions: translational motion along the vertical Z direction, rotational movement (lateral bending) and rotational movement ( flexion - extension). 

### Simulink simulations: 

#### Traction movement along Z 
&nbsp&nbsp&nbsp&nbsp&nbsp;                  Rotational Movement (Lateral Bending)
<p align="left">
  <img width="300" src="https://github.com/pranavs1911/TractionNeckBrace/blob/main/x2RPS1PRstract.gif">
  &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<img width="235" src="https://github.com/pranavs1911/TractionNeckBrace/blob/main/x2RPS1PRslatbending.gif">
  &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<img width="300" src="https://github.com/pranavs1911/TractionNeckBrace/blob/main/x2RPS1PRsflex.gif">
</p>

  

