% Initial Guess Values of Reaction Constants (sec ^ -1)

k1_not_plus = 2 * 10 ^ 2;
k1_not_minus = 2.6 * 10 ^ 10;
k2_plus = 2 * 10 ^ 11;
k2_minus = 1.3 * 10 ^ 4;
k3_plus = 1.0 * 10 ^ 12;
k3_minus = 5.8 * 10 ^ 2;
k4_plus = 1.5;
k4_minus = 6 * 10 ^ 5;
k5_plus = 8.0 * 10 ^ 4;
k5_minus = 1.0 * 10 ^ 11;
k6_plus = 6.0 * 10 ^ 7;
k6_minus = 2.4 * 10 ^ 8;


% Constants

T = 973;  % Absolute Temperature (K)
f_geo = 0.37;  % geometric factor
Ni_ss_density = 1.67 * 10^-9;  % Nickel Surface site Density (mol cm^-2)
A = 1;  % Electrode Surface Area (cm^-2)
z = 2;  % Number of electrons participating in the reaction
eta_var = 10; % Amplitude of input sinusoidal over-voltage (millivots)
eta_steady = 0;  % Steady state polarization over-voltage (millivolts)
F = 96500;  % Faraday's Constant (C mol^-1)
R = 8.3141; % Universal Gas Constant (J mol^-1 K^-1)
beta = 0.5; %  Symmetry Factor



% Overvoltage as a function of time

omega = 100; % hertz ----> omega = 10 or 1 gives error but 100 does not on solving the Thermodynamically consistent reaction kinetics model formulation
tau = ((2* pi)/(omega));
k1_plus = k1_not_plus * exp(beta * ((z*F)/(R*T)) * 0);
k1_minus = k1_not_minus * exp(-(1-beta) * ((z*F)/(R*T)) * 0);




% Functions used to solve the system of ode's


function d_theta_dt = odesystem(t2,theta,k1_plus,k1_minus,k2_plus,k2_minus,k3_plus,k3_minus,k4_plus,k4_minus,k5_plus,k5_minus,k6_plus,k6_minus)

  t2;
  theta_O = theta(1);
  theta_H = theta(2);
  theta_OH = theta(3);
  theta_H2O = theta(4);
  theta_Ni = 1 - (theta_O + theta_H + theta_OH + theta_H2O);


  d_theta_O_dt = (k1_plus * theta_Ni) - (k1_minus * theta_O) - (k2_plus * theta_H * theta_O) + (k2_minus * theta_OH * theta_Ni) - (k3_plus * theta_O * theta_H2O) + (k3_minus * theta_OH ^ 2);
  d_theta_H_dt = -(k2_plus * theta_O * theta_H) + (k2_minus * theta_OH * theta_Ni) - (k4_plus * theta_H * theta_OH) + (k4_minus * theta_H2O * theta_Ni) + (2 * k6_plus * theta_Ni ^ 2) - (2 * k6_minus * theta_H ^ 2);
  d_theta_OH_dt = (k2_plus * theta_O * theta_H) - (k2_minus * theta_OH * theta_Ni) + (2 * k3_plus * theta_O * theta_H2O) - (2* k3_minus * theta_OH ^ 2) - (k4_plus * theta_H * theta_OH) + (k4_minus * theta_H2O * theta_Ni);
  d_theta_H2O_dt = -(k3_plus * theta_O * theta_H2O) + (k3_minus * theta_OH ^ 2) + (k4_plus * theta_H * theta_OH) - (k4_minus * theta_H2O * theta_Ni) + (k5_plus * theta_Ni) - (k5_minus * theta_H2O); 

  d_theta_dt = [d_theta_O_dt ; d_theta_H_dt ; d_theta_OH_dt ; d_theta_H2O_dt];
  

end  


function d_theta_dt = odesystem2(t,theta2,R1,R2,R3,R4,R5,R6,cap_O,cap_H,cap_OH,cap_H2O,cap_Ni) %  here cap_i denotes the capacity value calculated using the equillibrium data

  F = 96500;  % Faraday's Constant (C mol^-1)
  R = 8.3141; % Universal Gas Constant (J mol^-1 K^-1)
  beta = 0.5; %  Symmetry Factor
  T = 973;
  z = 2;
  omega = 10;
  eta_var = 10; % Amplitude of input sinusoidal over-voltage (millivots)
  eta_steady = 0;  % Steady state polarization over-voltage (millivolts)
  eta = eta_steady + eta_var * sin (omega * t);


  e_f = exp(beta * ((z*F)/(R * T)) * eta);
  e_r = exp(-(1-beta) * ((z*F)/(R*T)) * eta);


  theta_O = theta2(1);
  theta_H = theta2(2);
  theta_OH = theta2(3);
  theta_H2O = theta2(4);
  theta_Ni = 1 - (theta_O + theta_H + theta_OH + theta_H2O);


  d_theta_O_dt = ((e_f/(R1 * cap_Ni))* theta_Ni) - ((e_r/(R1 * cap_O)) * theta_O) - ((1/(R2 * cap_O * cap_H)) * theta_H * theta_O) + ((1/(R2 * cap_OH * cap_Ni)) * theta_OH * theta_Ni) - ((1/(R3 * cap_O * cap_H2O)) * theta_O * theta_H2O) + ((1/(R3 * cap_OH ^ 2)) * theta_OH ^ 2);
  d_theta_H_dt = -((1/(R2 * cap_O * cap_H)) * theta_O * theta_H) + ((1/(R2 * cap_OH * cap_Ni)) * theta_OH * theta_Ni) - ((1/(R4 * cap_H * cap_OH)) * theta_H * theta_OH) + ((1/(R4 * cap_H2O * cap_Ni)) * theta_H2O * theta_Ni) + (2 * (1/(R6 * cap_Ni ^ 2)) * theta_Ni ^ 2) - (2 * (1/(R6 * cap_H ^ 2)) * theta_H ^ 2);
  d_theta_OH_dt = ((1/(R2 * cap_O * cap_H)) * theta_O * theta_H) - ((1/(R2 * cap_OH * cap_Ni)) * theta_OH * theta_Ni) + (2 * (1/(R3 * cap_O * cap_H2O)) * theta_O * theta_H2O) - (2 * (1/(R3 * cap_OH ^ 2)) * theta_OH ^ 2) - ((1/(R4 * cap_H * cap_OH)) * theta_H * theta_OH) + ((1/(R4 * cap_H2O * cap_Ni)) * theta_H2O * theta_Ni);
  d_theta_H2O_dt = -((1/(R3 * cap_O * cap_H2O)) * theta_O * theta_H2O) + ((1/(R3 * cap_OH ^ 2)) * theta_OH ^ 2) + ((1/(R4 * cap_H * cap_OH)) * theta_H * theta_OH) - ((1/(R4 * cap_H2O * cap_Ni)) * theta_H2O * theta_Ni) + ((1/(R5 * cap_Ni)) * theta_Ni) - ((1/(R5 * cap_H2O)) * theta_H2O); 

  d_theta_dt = [d_theta_O_dt ; d_theta_H_dt ; d_theta_OH_dt ; d_theta_H2O_dt];
  

end  






Wr = (k2_minus*k3_plus*k4_plus)/(k2_plus*k3_minus*k4_minus);  % Wegscheider Ratio


iterations = 0;


while abs(Wr - 1) > 10^-10 


    iterations = iterations + 1; 




    %  Step-1 : Solving the equillibrium model


    theta_initial = [0,0,0,0];
    tspan = [0 , 2*pi/omega];
    odeFunc = @(t2, theta) odesystem(t2, theta, k1_plus, k1_minus,k2_plus,k2_minus,k3_plus,k3_minus,k4_plus,k4_minus,k5_plus,k5_minus,k6_plus,k6_minus);
    [t2,theta] = ode15s(odeFunc,tspan,theta_initial);


    % Retrieve theta_i 's from the solution

    cap_O = theta(end, 1);
    cap_H = theta(end, 2);
    cap_OH = theta(end, 3);
    cap_H2O = theta(end, 4);
    cap_Ni = 1 - (cap_O + cap_H + cap_OH + cap_H2O);  % Compute theta_Ni based on the condition



    % Calculation of the resistances using the capacities calculated from the eqillibrium model 

    R1 = (1/k1_minus) * (cap_O ^ -1);
    R2 = (1/k2_plus) * (cap_O ^ -1) * (cap_H ^ -1);
    R3_guess = (1/k3_plus) * (cap_O ^ -1) * (cap_H2O ^ -1); % This is the initial guess value for R3
    R4_guess = (1/k4_plus) * (cap_OH ^ -1) * (cap_H ^ -1); % This is the initial guess value for R4
    R5 = (1/k5_minus) * (cap_H2O ^ -1);
    R6 = (1/k6_minus) * (cap_H ^ -2); 


    % using these resistances to solve the Thermodynamically consistent reaction kinetics model formulation


    theta_initial2 = [0,0,0,0];
    tspan2 = [0 , 2*pi/omega];
    odeFunc2 = @(t, theta2) odesystem2(t,theta2,R1,R2,R3,R4,R5,R6,cap_O,cap_H,cap_OH,cap_H2O,cap_Ni);
    [t,theta2] = ode15s(odeFunc2,tspan2,theta_initial2);
    
    
    % Retrieve theta_i 's from the solution
    theta_O = theta2(:, 1);
    theta_H = theta2(:, 2);
    theta_OH = theta2(:, 3);
    theta_H2O = theta2(:, 4);
    theta_Ni = 1 - (theta_O + theta_H + theta_OH + theta_H2O);  % Compute theta_Ni based on the condition


    % Plot the results


    figure;
    plot(t, theta_O, t, theta_H, t, theta_OH, t, theta_H2O, t, theta_Ni);
    legend('theta_O', 'theta_H', 'theta_OH', 'theta_H2O', 'theta_Ni');
    xlabel('Time t');
    ylabel('theta_i');
    title('Solution of Differential Equations with Condition theta_O + theta_H + theta_OH + theta_H2O + theta_Ni = 1');


    % PARAMETER ESTIMATION (R3 and R4) --->  Does not work 


    % Calculating the current:


    eta_var = 10; % Amplitude of input sinusoidal over-voltage (millivots)
    eta_steady = 0;  % Steady state polarization over-voltage (millivolts)
    eta = eta_steady + eta_var * sin (omega * t);

    k1_plus_t = k1_not_plus * exp(beta * ((z*F)/(R*T)) * eta);
    k1_minus_t = k1_not_minus * exp(-(1-beta) * ((z*F)/(R*T)) * eta);

    i_t = simplify((z * F * Ni_ss_density * A * f_geo) * ((k1_plus_t * theta_Ni)  -  (k1_minus_t * theta_O))); % i_t represents the current


    % Finding the Impeadence: Here we will use numerical integration using the "trapz" function

    sinfunc = sin(omega * t) ;
    cosfunc = cos(omega * t);

    %  Y represents the comples admittance

    Y_real = ((2 * A * f_geo)/(eta_var * tau)) * (trapz(t,(i_t * sinfunc))) ;
    Y_imaginary = ((2 * A * f_geo)/(eta_var * tau)) * (trapz(t,(i_t * cosfunc))) ;

    % Z represents the impeadence

    Y_square = simplify(Y_real ^ 2 + Y_imaginary ^ 2);
    Z_square = simplify(Y_square ^ -1);

    % Now we have an expression of impeadence in terms of R3 and R4. Now using the "fminsearch function", we will find the optimum value of R3 and R4.

    objective = @(R3,R4)Z_square;
    R_guess = [R3_guess,R4_guess];
    R_final = fminsearch(objective,R_guess);

    R3 = R_final(1);
    R4 = R_final(2);



  % Updating the K_i's 

    k1_minus = (1/R1) * (theta_O(end) ^ -1);
    k1_plus = (1/R1) * (theta_Ni(end) ^ -1);
    k2_plus = (1/R2) * (theta_O(end) ^ -1) * (theta_H(end)^ -1);
    k2_minus = (1/R2) * (theta_OH(end) ^ -1) * (theta_Ni(end) ^ -1);
    k3_plus = (1/R3) * (theta_O(end) ^ -1) * (theta_H2O(end)^ -1);
    k3_minus = (1/R3) * (theta_OH(end) ^ -2);
    k4_plus = (1/R4) * (theta_OH(end) ^ -1) * (theta_H(end) ^ -1);
    k4_minus = (1/R4) * (theta_H2O(end) ^ -1) * (theta_Ni(end) ^ -1);
    k5_plus = (1/R5) * (theta_Ni(end) ^ -1); % check once more i think there might be some mistake 
    k5_minus = (1/R5) * (theta_H2O(end) ^ -1);
    k6_plus = (1/R6) * (theta_Ni(end) ^ -2);
    k6_minus = (1/R6) * (theta_H(end) ^ -2);

    

     
    Wr = (k2_minus*k3_plus*k4_plus)/(k2_plus*k3_minus*k4_minus);  

end





    