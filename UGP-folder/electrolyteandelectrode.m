clear
close all
clc

% Complete cell code (electrode + eletrolyte)

L = 1500 * 10^-9;
M = 320 * 10^-9;
cellsize = L + M;
nx = 183;
nt = 3600*80 + 1;
dom_size = 1820* 10^-9;     % domsize = cellsize
timedom_size = 3600;
dx = dom_size/(nx - 1);
dt = timedom_size/(nt-1);
D_Liplus = 0.9 * 10^-15;
D_nminus = 5.1 * 10^-15;
D_final  = (2 * D_Liplus * D_nminus) / (D_Liplus + D_nminus);
D_Li = 1.76 * 10^-15;
delta = 0.18;
a_not = 6.01 * 10^4;
a_not_licoo2 = 2.33 * 10^4;
k_r = 0.90 * 10^-8;
k_d = ((k_r * a_not * delta * delta) / (1 - delta)) ;
A = 10^-4;
F = 96500;
I = 60 * 10^-6;


alpha1 = (D_final * (dt/(dx*dx))); 
beta1 = (I/(2 * F * A * D_Liplus));
alpha2 = (D_Li * (dt/(dx*dx)));
beta2 = (I/(F * A * D_Li));
alpha3 = ((2*(D_Liplus))/(D_Li));
a = zeros(nt,nx);



a (1,1:151) = delta * a_not;
a (1,152:183) = a_not_licoo2;



% Initialize figure
legend;
xlabel('Position');
ylabel('Time');
title('Concentration Profile for complete Cell (Electrode + Electrode)');
hold on;



for i = 2 : nt

    for j = 2 : 150

       a(i,j) = a(i-1,j) + alpha1 * (a(i-1,j+1) - 2*a(i-1,j)+ a(i-1,j-1)) + dt * ((k_d * a_not) - (k_d * a(i-1,j)) - (k_r * a(i-1,j) * a(i-1,j)));

    end  


    for j = 152 : nx-1

        a(i,j) = a(i-1,j) + alpha2 * (a(i-1,j+1) - 2*a(i-1,j)+ a(i-1,j-1));

    end


    a(i,1) = a(i,2) - beta1 * dx;   % left boundary

    a(i,nx) = a(i,nx-1);    % right boundary

    a(i,151) = ((a(i,152) + alpha3 * a(i,150))/(1 + alpha3));    % common boundary between electrolyte and electrode



    if  i == 1000 || i == 5000 || i == 10000|| i == 25000 || i == 50000 || i == 100000 || i == 150000 || i == 200000 || i == 250000 ||i == 288000
        
        
        plot((10^-8)*(1:nx), a(i, :), 'DisplayName', ['Time Step ', num2str(i)]);


    end

end    


hold off;
legend;
xlim((10^-8)*[1, nx]);
xlabel('Position');
ylabel('Concentration');
title('Concentration Profile for complete Cell (Electrode + Electrode)');










