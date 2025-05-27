% simulate.m
clear; clc; close all

%parameters
tspan = [0 20]; 

%initial angles
phib0  = deg2rad(0);
theta0 = deg2rad(10);

q_ind0 = [phib0; theta0];

%initial angular velocities
phib_dot0 = 1;
theta_dot0 = 0;
qdot_ind0 = [phib_dot0; theta_dot0];

%State vector: x = [phib; theta; phib_dot; theta_dot] 
x0 = [q_ind0; qdot_ind0];

options = odeset('RelTol',1e-6, 'AbsTol',1e-8);

[T, X] = ode45(@dynamics, tspan, x0, options);

%trajectory
phib_dot = X(:,3);     
theta    = X(:,2);      

v = 0.05 * phib_dot; 

xPos = cumtrapz(T, v .* cos(theta));
yPos = cumtrapz(T, v .* sin(theta));


margin = 0.1;  
xmin = min(xPos) - margin;
xmax = max(xPos) + margin;
ymin = min(yPos) - margin;
ymax = max(yPos) + margin;

r  = 0.05;    
L1 = 0.28;    
L2 = 0.28;   

phib_dot  = X(:,3);  
theta     = X(:,2);  
theta_dot = X(:,4);  
psidot = -( L2*theta_dot + phib_dot.*r.*sin(theta) ) ...
         ./ ( L2 + L1*cos(theta) );


psi = cumtrapz(T, psidot);


v = r * phib_dot;


xPos = cumtrapz(T, v .* cos(psi));
yPos = cumtrapz(T, v .* sin(psi));

%plot
figure;
plot(xPos, yPos, 'LineWidth',1.5);
axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('Robot Trajectory');

figure;
%phi
subplot(2,1,1);
plot(T, rad2deg(X(:,1)), 'LineWidth',1.5);
xlabel('Time (s)');
ylabel('φ (deg)');
title('φ vs. Time');
grid on;
%phi velocity
subplot(2,1,2);
plot(T, X(:,3), 'LineWidth',1.5);
xlabel('Time (s)');
ylabel('φ̇ (rad/s)');
title('φ̇ vs. Time');
grid on;


figure;
%theta 
subplot(2,1,1);
plot(T, rad2deg(X(:,2)), 'LineWidth',1.5);
xlabel('Time (s)');
ylabel('θ (deg)');
title('θ vs. Time');
grid on;
%theta velocity
subplot(2,1,2);
plot(T, X(:,4), 'LineWidth',1.5);
xlabel('Time (s)');
ylabel('θ̇ (rad/s)');
title('θ̇ vs. Time');
grid on;

%energy
E = zeros(length(T),1);
for i = 1:length(T)
    theta_i     = X(i,2);
    phib_dot_i  = X(i,3);
    theta_dot_i = X(i,4);

    M11 = (4917*cos(theta_i) + 4932) / (9800*(cos(theta_i)+1));
    M12 = (4932*sin(theta_i) + 4917*cos(theta_i)*sin(theta_i)) ...
          / (3500*(cos(theta_i)+1)^2);
    M22 = -(4917*cos(theta_i)^2 - 4932) / (625*(cos(theta_i)+1)^2);
    M_i = [ M11, M12;
            M12, M22 ];

    qdot = [phib_dot_i; theta_dot_i];
    E(i) = 0.5 * qdot' * M_i * qdot;
end

% plot
figure
plot(T,E,'LineWidth',1.5);     
grid on
xlabel('Time (s)')
ylabel('Kinetic Energy (J)')
title('System Energy vs. Time')
xlim([0 10])
xticks(0:1:10)


