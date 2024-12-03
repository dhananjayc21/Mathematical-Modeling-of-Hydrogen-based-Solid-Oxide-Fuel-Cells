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


  %  Step-1 : Solving the equillibrium model


    theta_initial = [0,0,0,0];
    tspan = [0 , 2*pi/omega];
    odeFunc = @(t, theta) odesystem(t, theta, k1_plus, k1_minus,k2_plus,k2_minus,k3_plus,k3_minus,k4_plus,k4_minus,k5_plus,k5_minus,k6_plus,k6_minus);
    [t,theta] = ode15s(odeFunc,tspan,theta_initial);


    % Retrieve theta_i 's from the solution

    theta_O = theta(:, 1);
    theta_H = theta(:, 2);
    theta_OH = theta(:, 3);
    theta_H2O = theta(:, 4);
    theta_Ni = 1 - (cap_O + cap_H + cap_OH + cap_H2O);  % Compute theta_Ni based on the condition



    figure;
    plot(t, theta_O, t, theta_H, t, theta_OH, t, theta_H2O, t, theta_Ni);
    legend('theta_O', 'theta_H', 'theta_OH', 'theta_H2O', 'theta_Ni');
    xlabel('Time t');
    ylabel('theta_i');
    title('Solution of Differential Equations with Condition theta_O + theta_H + theta_OH + theta_H2O + theta_Ni = 1');



