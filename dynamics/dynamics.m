function xdot = dynamics(t, x)

    syms psib theta phib real
    syms psib_dot theta_dot phib_dot real

    % x = [phib; theta; phib_dot; theta_dot]

    q = x(1:2);
    qdot = x(3:4);
    fullState = {q(1), q(2), qdot(1), qdot(2)};
    
    for k = 1:length(fullState)
        fullState{k} = double(fullState{k});
    end
   
    
    M_reduced = [                        (4917*cos(theta) + 4932)/(9800*(cos(theta) + 1)), (4932*sin(theta) + 4917*cos(theta)*sin(theta))/(3500*(cos(theta) + 1)^2);
                (4932*sin(theta) + 4917*cos(theta)*sin(theta))/(3500*(cos(theta) + 1)^2),                     -(4917*cos(theta)^2 - 4932)/(625*(cos(theta) + 1)^2)];


    C_reduced = [(315*phib_dot - 315*phib_dot*cos(theta)^2 - 315*phib_dot*cos(theta)^3 + ...
        315*phib_dot*cos(theta) + 1884*theta_dot*sin(theta) + ...
        1764*theta_dot*cos(theta)*sin(theta))/(156800*(cos(theta) + 1)^2), (3*theta_dot*(3273*cos(theta) + 3293))/(7000*(cos(theta) + 1)^2);
        (3*(608*theta_dot - 588*theta_dot*cos(theta)^2 + 105*phib_dot*sin(theta)^3 - ...
        20*theta_dot*cos(theta)))/(28000*(cos(theta) + 1)^2), (3*theta_dot*sin(theta)*(1639*cos(theta) + 1644))/(625*(cos(theta) + 1)^3)];

    allVars = [phib, theta, phib_dot, theta_dot];
    M_func = matlabFunction(M_reduced, 'Vars', allVars);
    C_func = matlabFunction(C_reduced, 'Vars', allVars);
 
    M_val = M_func(fullState{:});
    C_val = C_func(fullState{:});

    %commanded torques
    tau_phib_cmd  = 0;      
    tau_theta_cmd = 0;      

    %damping coefficients
    b_phi   = 0;    % 0.05 wheel damping
    b_theta = 0;    % 0.2 joint damping 

    phi_dot   = x(3);    
    theta_dot = x(4);    

    % Total torques with damping
    tau_phib  = tau_phib_cmd  - b_phi   * phi_dot;
    tau_theta = tau_theta_cmd - b_theta * theta_dot;

    F_val = [tau_phib; tau_theta];

    qdd = M_val \ (F_val - C_val * qdot);
    
    xdot = [qdot; qdd];
end
