syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot;
syms m1 m2 l1 l2 r1 r2 I1 I2 g T1 T2;
syms theta_1_t_0 theta_1_tf theta_2_t_0 theta_2_tf theta_1_dt_0 theta_1_dtf theta_2_dt_0 theta_2_dtf
syms j1_a0 j1_a1 j1_a2 j1_a3 j2_a0 j2_a1 j2_a2 j2_a3;
syms t0 tf t



K = 0.5*m1*(theta1_dot*r1)^2 + 0.5*I1*theta1_dot^2 + 0.5*m2*((theta1_dot*l1)^2+(r2*(theta1_dot+theta2_dot))^2 + 2*l1*theta1_dot*r2*(theta1_dot+theta2_dot)*cos(theta2))+0.5*I2*(theta1_dot+theta2_dot)^2;

P = m1*g*r1*cos(theta1)+m1*g*l1*cos(theta1)+m2*g*r2*cos(theta1+theta2);
L = K - P;
d_dt_l = jacobian(L,[theta1,theta2,theta1_dot,theta2_dot]);
d_dt_l_1 = jacobian(d_dt_l(3),[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];
d_dt_l_2 = jacobian(d_dt_l(4),[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];

eq1 = d_dt_l_1-d_dt_l(1)-T1;
eq2 =d_dt_l_2-d_dt_l(2)-T2;

solution = solve([eq1==0 , eq2==0] ,[theta1_ddot,theta2_ddot]);
solution.theta1_ddot
solution.theta2_ddot

g1 = simplify(subs((eq1), [theta1_dot, theta2_dot, theta1_ddot, theta2_ddot],[0,0,0,0])) + T1
M1 = simplify(subs((eq1), [theta1_dot, theta2_dot], [0,0]) - g1)
C1 = simplify(subs(eq1) - M1 - g1)

g2 = simplify(subs((eq2), [theta1_dot, theta2_dot, theta1_ddot, theta2_ddot],[0,0,0,0])) + T2
M2 = simplify(subs((eq2), [theta1_dot, theta2_dot], [0,0]) - g2)
C2 = simplify(subs(eq2) - M2 - g2)

M = jacobian([M1;M2], [theta1_ddot; theta2_ddot])
g = [g1;g2];
C = [0 C1; C2 0];


tf = subs(tf,10);
t0 = subs(t0, 0);

theta_1_t_0 = subs(theta_1_t_0, pi);
theta_2_t_0 = subs(theta_2_t_0, pi/2);

theta_1_tf = subs(theta_1_tf, 0);
theta_2_tf = subs(theta_2_tf, 0);

theta_1_dt_0 = subs(theta_1_dt_0, 0);
theta_2_dt_0 = subs(theta_2_dt_0, 0);

theta_1_dtf = subs(theta_1_dtf, 0);
theta_2_dtf = subs(theta_2_dtf, 0);

j1_e1 = j1_a0 + (j1_a1 * t0) + (j1_a2 * (t0)^2) + (j1_a3 * (t0)^3) - theta_1_t_0;
j1_e2 = j1_a1 + (2 * (j1_a2 * t0)) + (3 * (j1_a3 * (t0)^2)) - theta_1_dt_0;
j1_e3 = j1_a0 + (j1_a1 * tf) + (j1_a2 * (tf)^2) + (j1_a3 * (tf)^3) - theta_1_tf;
j1_e4 = j1_a1 + (2 * (j1_a2 * t0)) + (3 * (j1_a3 * (t0)^2)) - theta_1_dtf;


j2_e1 = j2_a0 + (j2_a1 * t0) + (j2_a2 * (t0)^2) + (j2_a3 * (t0)^3) - theta_2_t_0;
j2_e2 = j2_a1 + (2 * (j2_a2 * t0)) + (3 * (j2_a3 * (t0)^2)) - theta_2_dt_0;
j2_e3 = j2_a0 + (j2_a1 * tf) + (j2_a2 * (tf)^2) + (j2_a3 * (tf)^3) - theta_2_tf;
j2_e4 = j2_a1 + (2 * (j2_a2 * t0)) + (3 * (j2_a3 * (t0)^2)) - theta_2_dtf;

solve([j1_e1==0, j1_e2==0, j1_e3==0, j1_e4==0], [j1_a0, j1_a1, j1_a2, j1_a3]);
solve([j2_e1==0, j2_e2==0, j2_e3==0, j2_e4==0], [j2_a0, j2_a1, j2_a2, j2_a3]);



A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];

lambda = [ -0.6 -0.08 -6 -10]; 
K = place(A,B,lambda);

jq = [theta1; theta1_dot; theta2; theta2_dot];

jq_dot = [j1_a0 + j1_a1 * t + j1_a2 * t^2 + j1_a3 * t^3;
          j1_a1 + 2 * (j1_a2 * t) + 3 * (j1_a3 * (t)^2);
          j2_a0 + j2_a1 * t + j2_a2 * t^2 + j2_a3 * t^3;
          j2_a1 + 2 * (j2_a2 * t) + 3 * (j2_a3 * (t)^2)]

jq_ddot = [2 * j1_a2 + 6 * j1_a3 * t;
           2 * j2_a2 + 6 * j2_a3 * t]

Vir_U = M * (-K *(jq - jq_dot) + jq_ddot) + C + g