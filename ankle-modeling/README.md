# Ankle Modeling

Model the ankle stiffness analyticaly.
Assumes the ankle is a spherical joint with many linear springs and dampers in parallel.
The torque and impedance around the ankle can be computed analytically given:
 - Pose of the spring/dampers
 - Natural length of spring
 - Stiffness and damping coefficients of spring/dampers
 - Ankle angle

## Result

 - The impedance is nonlinear with the ankle angle.
 - The local ankle stiffness and damping takes the form of symmetric matrices