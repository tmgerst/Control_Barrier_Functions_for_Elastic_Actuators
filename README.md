# Control Barrier Functions for Learning-Based Regulation of Elastic Actuators with State Constraints

The code contained in this repo constitutes the first part of my bachelor thesis at the chair of Information-Oriented Control at TUM. 
It performs simulations on a two degree of freedom planar robot subject to constraints. 

Motivated by the huge importance of patient safety in medical applications, constraints in a robot's workspace or jointspace have to be ensured, since constraint violation in these scenarios can damage both robot and human.
One way to formally guarantee constraint satisfaction in controlled dynamical systems is the use of control barrier functions (CBFs). 
These functions create a sort of safety filter on the systems's control inputs: only control inputs are allowed that result in the system satisfying its constraints. 
And if a control input would result in constraint violation, the filter modifies it to ensure constraint satisfaction. 

## How to run
Simply run `main`.

It includes a rather straightforward simulation, where the parameters can be chosen right at the top. You can also choose from a handful of 
output options, including a small animation. Controller parameters can be changed in the appropriate functions in the control/` directory if needed, but the default should work fine.

## Code structure
| Directory / File | Description | 
| ---------------- | ----------- |
| control | A couple of controllers for the robot, which were tried during the project. |
| models | Some robot models with different properties, as well as a function to compute the inverse kinematics for one of them. |
| utilities | Utility functions for copying large amounts of variables around or setting initial states. |
| visualization | Two functions for drawing and animating outputs of the simulation. |
| `enforce_constraints.m` | Optimization for the control inputs based on CBFs encoding the constraints set in `main`. | 
| `main.m` | Main simulation and entry point into the project. |

## Requirements
A standard MATLAB installation without any toolboxes should work fine. 

## Project Description 
For this project, the focus was on the arm of an exoskeleton used in rehabilitation. The arm is driven by an elastic actuator (EA), where the driving motor
is decoupled from its link by a spring and damping elements. Hence, a motor motion is not directly transfered to the link, posing a challenge for its 
control. Since such actuators feature elastic elements however, systems driven by EAs are inherently compliant to patients' movements, making them a 
good choice in medical applications. 
The exoskeleton has to adhere to link angle constraints imposed by the patient, i.e the elbow angle might not exceed 60°. This has to even be ensured
when the reference trajectory violates this constraint. That's where the CBF control comes in, which modifies the nominal trajectory tracking control
command for the motor, such that the link angle always stays within its constraint, ensuring that the patient is not harmed by the robot. 



## Further work
In addition to the simulation in this repo, there also exists real-time LabView code for experiments on a real exoskeleton at the chair of Information-Oriented Control, which I have also written as part of this project. Results are presented in my thesis, shoot me a message in case you are interested!
