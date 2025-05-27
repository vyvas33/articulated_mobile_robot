    clear; clc;
    
    %symbolic variables
    
    syms xb yb psib theta phib phit real
    syms xb_dot yb_dot psib_dot theta_dot phib_dot phit_dot real
    
    L1 = 0.28;   % Distance from joint to front module wheels
    L2 = 0.28;   % Distance from joint to trailer module wheels   
    r  = 0.05;   % wheel radius
    
    mb = 4.5;    % mass of base module (front)
    mt = 4.5;    % mass of trailer module
    mw = 0.12;   % mass of wheels 
    Ib = 0.012;  % moment of inertia of base module (front)
    It = 0.012;  % moment of inertia of trailer module      
    Iw = 0.01;   % moment of inertia of wheels
    
    
    y   = [xb;     yb;     psib;     theta;     phib;     phit];
    y_d = [xb_dot; yb_dot; psib_dot; theta_dot; phib_dot; phit_dot];
    
    % kinetic energy
    
    T = (mt*((L2*cos(psib + theta)*(psib_dot + theta_dot) - yb_dot + L1*psib_dot*cos(psib))^2 + ...
        (xb_dot + L2*sin(psib + theta)*(psib_dot + theta_dot) + L1*psib_dot*sin(psib))^2))/2 + ...
        (It*(psib_dot + theta_dot)^2)/2 + ...
        (Ib*psib_dot^2)/2 + ...
        mw*(phib_dot^2 + phit_dot^2) + ...
        (mb*(xb_dot^2 + yb_dot^2))/2;
    
    %   M(i,j) = d^2 T / d(qdot_i) d(qdot_j)
    
    n = length(y); 
    M_sym = sym('M', [n n]);
    
    for i = 1:n
        for j = 1:n
            M_sym(i,j) = diff(diff(T, y_d(i)), y_d(j));
        end
    end
    
    M_sym = simplify(M_sym);
    
    %   Compute the Coriolis/Centrifugal Matrix C
    %   Using the Christoffel symbol approach:
    %   C_{i,j} * qdot_j = sum_{k} [C_{i,j,k} * qdot_j * qdot_k],
    %   where
    %   C_{i,j,k} = 0.5 * ( dM(i,j)/dq(k) + dM(i,k)/dq(j) - dM(j,k)/dq(i)).
    
    C_sym = sym('C', [n n]);
    
    for i = 1:n
        for j = 1:n
            c_temp = sym(0);
            for k = 1:n
                c_ijk = 0.5*( diff(M_sym(i,j), y(k)) ...
                            + diff(M_sym(i,k), y(j)) ...
                            - diff(M_sym(j,k), y(i)) );
                c_temp = c_temp + c_ijk*y_d(k);
            end
            C_sym(i,j) = c_temp;
        end
    end
    
    C_sym = simplify(C_sym);
     
    % Inorder to incorporate the kinematic constraints, we incorporate lagrange
    % multipliers into the lagrange equations A(q)'*lambda = 0  
    
    J =  [                        r*cos(psib),                                           0;
                                  r*sin(psib),                                           0;
         -(r*sin(theta))/(L2 + L1*cos(theta)),                    -L2/(L2 + L1*cos(theta));
                                            0,                                           1;
                                            1,                                           0;
    (L1 + L2*cos(theta))/(L2 + L1*cos(theta)), (L1*L2*sin(theta))/(r*(L2 + L1*cos(theta)))];
    
     
    J_dot = [                                                          0,                                                                 0
                                                                      0,                                                                 0
             -(r*theta_dot*(L1 + L2*cos(theta)))/(L2 + L1*cos(theta))^2,              -(L1*L2*theta_dot*sin(theta))/(L2 + L1*cos(theta))^2
                                                                      0,                                                                 0
                                                                      0,                                                                 0
            (theta_dot*sin(theta)*(L1^2 - L2^2))/(L2 + L1*cos(theta))^2, (L1*L2*theta_dot*(L1 + L2*cos(theta)))/(r*(L2 + L1*cos(theta))^2)];
    
    
    M_reduced = J'*M_sym*J;
    disp(simplify(M_reduced));
    
    C_reduced = J'*M_sym*J_dot + J'*C_sym*J;
    C_reduced = subs(C_reduced,psib_dot,-(L2*theta_dot + phib_dot*r*sin(theta))/(L2 + L1*cos(theta)));
    disp(simplify(C_reduced));
