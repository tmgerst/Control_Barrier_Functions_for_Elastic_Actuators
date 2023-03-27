# Control Barrier Functions for Learning-Based Regulation of Elastic Actuators with State Constraints

## Bachelor Thesis at the Chair of Information-Oriented Control of the Technical University of Munich

The code contained in this repo performs simulations on a two degree of freedom planar robot subject to constraints. 
Motivated by the thought that in many medical applications patient safety is of the utmost importance, 
constraints in a robot's workspace or jointspace have to be ensured, since constraint violation in these scenarios can damage both robot and human.
One way to formally guarantee constraint satisfaction in controlled dynamical systems is the use of control barrier functions (CBFs), which create a safety 
filter on the systems's control inputs, allowing only the inputs that ensure that the system remains within its constraints.

## Project Description 
For this project, the focus lies on the arm of an exoskeleton used in rehabilitation. The arm is driven by an elastic actuator (EA), where the driving motor
is decoupled from its link by a spring and damping elements. Hence, a motor motion is not directly transfered to the link, posing a challenge from a 
control perspective. Since they feature elastic elements however, systems driven by EAs are inherently compliant to patients' movements, making them a 
good choice in medical applications. 
The exoskeleton has to adhere to link angle constraints imposed by the patient, i.e the elbow angle might not exceed 60°. This has to even be ensured
when the reference trajectory violates this constraint. That's where the CBF control comes in, which modifies the nominal trajectory tracking control
command for the motor, such that the link angle always stays within its constraint, ensuring that the patient is not harmed by the robot. 

For more theory as well as a small animation, check out the presentation. The animation can be found on page 7 and is activated with a click.

## Code structure
The main file includes a rather straightforward simulation, where the parameters can be chosen right at the top. You can also choose from a handful of 
output options, including a small animation. Controller parameters can be changed in the appropriate functions in the 'control' directory if needed, but the default should work fine.

## Further work
The simulations have been  performed to tune the CBFs to the presence of dampings inside the joints. 
Afterwards, the CBF control has also been implemented in LabView in order to flash the exoskeleton at the chair of Information-Oriented Control and
experiments were performed (see the presentation for some results). 
