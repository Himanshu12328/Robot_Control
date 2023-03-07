syms t t_0 tf a_0 a_1 a_2 a_3 theta_t_0 theta_d_t_0 theta_tf theta_d_tf 'real'


A = [ 1, t_0, t_0^02, t_0^3;
      0, 1 , 2*t_0, 3*t_0^2;
      1, tf, tf^2, tf^3;
      0, 1, 2*tf, 3*tf^2];

a = [a_0; 
     a_1; 
     a_2; 
     a_3];

b = [theta_t_0;
     theta_d_t_0;
     theta_tf;
     theta_d_tf];
A1 = [1, 0, 0, 0;
      0, 1, 0, 0;
      1, 10, 100, 1000;
      0, 1, 20, 300];

b1 = [pi;
       0;
       0;
       0];

a1 = pinv(A1) * b1

theta_1_t = a1(1) + a1(2)*t + a1(3)*t^2 + a1(4)*t^3
theta_1_d_t = a1(2) + 2*a1(3)*t + 3*a1(4)*t^2

A2 = [1, 0, 0, 0;
      0, 1, 0, 0;
      1, 10, 100, 1000;
      0, 1, 20, 300];
b2 = [pi/2;
         0;
         0;
         0];
a2 = pinv(A2) * b2

theta_2_t = a2(1) + a2(2)*t + a2(3)*t^2 + a2(4)*t^3
theta_2_d_t = a2(2) + 2*a2(3)*t + 3*a2(4)*t^2

theta_1_dd_t = 2*a1(3) + 6*a1(4)*t
theta_2_dd_t = 2*a2(3) + 6*a2(4)*t

