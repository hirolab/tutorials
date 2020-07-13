This document explains how to use matlab's linear solver to fit parameters.

The example code uses a made up spring mass damper who is given a random trajectory (displacement) in a given time vector. 

Then, the displacement derivatives are calculated (velocity and displacement).

Values for stiffness(K), damping(B), inertia(I), hyesteresis(C) and force offset(G) are all defined. We expect our parameter predictor to return those values.

The equation of motion is developed and solved for force, the only parameter left without data.



Scenario: In normal experiments we can measure displacement (using the mocap or IMU) and we can measure force (usign a load cell or the force plate).

Our objective is the characterize the stiffness, damping, inertia, and hysteresis, as well as any force offset so it doesn't affect our data.

Since we know how these parameters are related to each other, we can develop the equation of motion:



F = x*K + dx*B + ddx*I + C*sign(dx) + G



Now we can use matrices to solve the parameters. The matrix arrangement is as follows:

The A matrix is made up of concatanated column vectors of data over time:



A = [x(1) dx(1) ddx(1) sign(dx(1) 1

     x(2) dx(2) ddx(2) sign(dx(2) 1

     ...

     x(n) dx(n) ddx(n) sign(dx(n) 1];

    

Also, F must be presented in column form!



Now Q is the column vector of the parameter values we want to obtain, conveniently named by the parameters we want to find.

Normally, mathematically we can solve a system of equation in the form of:



[F(1)   [x(1) dx(1) ddx(1) sign(dx(1) 1     [Q(1)  

 F(2) =  x(2) dx(2) ddx(2) sign(dx(2) 1   *  Q(2)   

 ...     ...                                        

 F(n)]   x(n) dx(n) ddx(n) sign(dx(n) 1];    Q(n)];

 

Luckily, matlab has such function that by writing Q = A\B, it returns a column vector of all the parameters solved usign least squares! 

Thus we successfully identified parameters. Remember that the sample code is used to validate the system, but it should be applicable to ankle and prosthesis parameter ID

