clear
close all
clc  

% Electrolyte only code

nx = 21;
nt = 3601;
dom_size = 1820 * 10^-9;
timedom_size = 3600;
dx = (dom_size/(nx - 1));
dt = (timedom_size/(nt - 1));
D_Liplus = 0.9 * 10^-15;
D_nminus = 5.1 * 10^-15;
D_final  = (2 * D_Liplus * D_nminus) / (D_Liplus + D_nminus);
delta = 0.18;
a_not = 6.01 * 10^4;
k_r = 0.90 * 10^-8;

k_d = ((k_r * a_not * delta * delta) / (1 - delta)) ;
A = 10^-4;
F = 96500;
I = 60 * 10^-6;

alpha = (D_final * (dt/(dx*dx))); 
beta = (I/(2 * F * A * D_Liplus));

a = zeros(nt,nx);

a (1,:) = delta * a_not;

% Initialize figure
legend;
xlabel('Position');
ylabel('Time');
title('Concentration Profile for electrolyte only');
hold on;

for i = 2 : nt


    for j = 2 : nx-1


       a(i,j) = a(i-1,j) + alpha * (a(i-1,j+1) - 2*a(i-1,j)+ a(i-1,j-1)) + dt * ((k_d * a_not) - (k_d * a(i-1,j)) - (k_r * a(i-1,j) * a(i-1,j)));

    end  


    a(i,nx) = a(i,nx-1) - beta * dx;

    a(i,1) = a(i,2) + beta * dx ;


     if  i == 100 || i == 200 || i == 300|| i == 400 || i == 500 || i == 600 || i == 700 || i == 800 || i == 900 ||i == 1000
        
        
        plot((10^-8)*(1:nx), a(i, :), 'DisplayName', ['Time Step ', num2str(i)]);
        

    end

end    


hold off;
legend;
xlim((10^-8)*[1, nx]);
xlabel('Position');
ylabel('Concentration');
title('Concentration Profiles for electrolyte only');


% Plotting
figure;
[X, T] = meshgrid(linspace(0, dom_size, nx), linspace(0, timedom_size, nt));
surf(T, X, a);









