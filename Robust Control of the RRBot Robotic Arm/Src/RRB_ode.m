%ODE function
function dX = RRB_ode(t,X)

    m1 = 1;
    m2 = 1;
    l1 = 1;
    %l2 = 1;
    r1 = 0.45;
    r2 = 0.45;
    I1 = 0.084;
    I2 = 0.084;
    g = 9.81;

    dX = zeros(4,1);
    X = num2cell(X);
    [theta1, theta2, theta1dot, theta2dot] = deal(X{:});
    
%     K = [26.2200   -7.6202    6.7730   -4.8545;
%          6.8728    0.9502    2.2068   -1.2374];
%         
%     T1 = -K(1,:) * [theta1 theta2 theta1dot theta2dot]';
%     T2 = -K(2,:) * [theta1 theta2 theta1dot theta2dot]';

    j1_a0 = pi;
    j1_a1 = 0;
    j1_a2 = -(3*pi)/100;
    j1_a3 = pi/500;
    
    j2_a0 = pi/2;
    j2_a1 = 0;
    j2_a2 = -(3*pi)/200;
    j2_a3 = pi/1000;

    A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
    B = [0 0; 0 0; 1 0; 0 1];
    
    lambda = [-0.6 -0.08 -6 -10]; 
    K = place(A,B,lambda);

    q = [ theta1; theta2; theta1dot; theta2dot];

    q_d = [j1_a0 + j1_a1 * t + j1_a2 * t^2 + j1_a3 * t^3;
            j2_a0 + j2_a1 * t + j2_a2 * t^2 + j2_a3 * t^3;
            j1_a1 + 2 * j1_a2 * t + 3 * j1_a3 * t^2;
            j2_a1 + 2 * j2_a2 * t + 3 * j2_a3 * t^2];

    q_dd = [  2 * j1_a2 + 6 * j1_a3 * t;
                2 * j2_a2 + 6 * j2_a3 * t;];

    G = [- g*m2*r2*sin(theta1 + theta2) - g*l1*m1*sin(theta1) - g*m1*r1*sin(theta1);
         -g*m2*r2*sin(theta1 + theta2)];

    C = [-l1*m2*r2*theta2dot*sin(theta2)*(2*theta1dot + theta2dot);
         l1*m2*r2*theta1dot^2*sin(theta2)];
    M = [m1*r1^2 + I1 + I2 + (m2*(2*l1^2 + 4*cos(theta2)*l1*r2 + 2*r2^2))/2, I2 + (m2*(2*r2^2 + 2*l1*cos(theta2)*r2))/2;
         I2 + (m2*(2*r2^2 + 2*l1*cos(theta2)*r2))/2,                               m2*r2^2 + I2];
    
    Vir_U = M*(-K * (q - q_d) + q_dd ) + C + G;

    T1 = Vir_U(1);
    T2 = Vir_U(2);


    dX(1) = theta1dot;
    dX(2) = theta2dot;
    dX(3) = (I2*T1 - I2*T2 + T1*m2*r2^2 - T2*m2*r2^2 + l1*m2^2*r2^3*theta1dot^2*sin(theta2) + l1*m2^2*r2^3*theta2dot^2*sin(theta2) + I2*l1*g*m1*sin(theta1) - l1*T2*m2*r2*cos(theta2) + I2*g*m1*r1*sin(theta1) + 2*l1*m2^2*r2^3*theta1dot*theta2dot*sin(theta2) + l1^2*m2^2*r2^2*theta1dot^2*cos(theta2)*sin(theta2) - l1*g*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1dot^2*sin(theta2) + I2*l1*m2*r2*theta2dot^2*sin(theta2) + l1*g*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1dot*theta2dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
    dX(4) = -(I2*T1 - I1*T2 - I2*T2 - l1^2*T2*m2 - T2*m1*r1^2 + T1*m2*r2^2 - T2*m2*r2^2 + l1*m2^2*r2^3*theta1dot^2*sin(theta2) + l1^3*m2^2*r2*theta1dot^2*sin(theta2) + l1*m2^2*r2^3*theta2dot^2*sin(theta2) - l1^2*g*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*l1*g*m1*sin(theta1) + l1*T1*m2*r2*cos(theta2) - 2*l1*T2*m2*r2*cos(theta2) + I2*g*m1*r1*sin(theta1) + 2*l1*m2^2*r2^3*theta1dot*theta2dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2dot^2*cos(theta2)*sin(theta2) - l1*g*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I1*l1*m2*r2*theta1dot^2*sin(theta2) + I2*l1*m2*r2*theta1dot^2*sin(theta2) + I2*l1*m2*r2*theta2dot^2*sin(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + l1*g*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1dot*theta2dot*cos(theta2)*sin(theta2) + l1^2*g*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1dot*theta2dot*sin(theta2) + l1*g*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
end
