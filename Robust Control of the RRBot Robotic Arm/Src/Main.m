syms m1_hat m2_hat I1_hat I2_hat 'real'

j1_a0 = pi;
j1_a1 = 0;
j1_a2 = -(3*pi)/100;
j1_a3 = pi/500;

j2_a0 = pi/2;
j2_a1 = 0;
j2_a2 = -(3*pi)/200;
j2_a3 = pi/1000;

m1_hat = 0.75;
m2_hat = 0.75;

I1_hat = 0.063;
I2_hat = 0.063;

A = [0, 0, 1, 0;
     0, 0, 0, 1;
     0, 0, 0, 0;
     0, 0, 0, 0]
B = [0, 0;
     0, 0;
     1, 0;
     0, 1]

lambda = [-3, -3, -4, -4]

K = place(A, B, lambda)
kp = K(:,1:2)
kd = K(:,3:4)

A_cl = [0, 0, 1, 0;
        0, 0, 0, 1;
      -12, 0, -7, 0;
        0, -12, 0, -7]
Q = eye(4)
P = lyap(A_cl' , Q)
rho = 1;
phi = 0.05;

tf = 10;
xx0 = [deg2rad(200), deg2rad(125), 0, 0];
[tf, X] = ode45(@(t, x) RRB_ode(t, x, K, P, rho, phi), [0, tf], xx0);

U = [];

for i = 1:length(tf)
    xx1(i:1) = j1_a0 + j1_a1*tf(i) + j1_a2*tf(i)^2 + j1_a3*tf(i)^3;
    xx2(i:1) = j2_a0 + j2_a1*tf(i) + j2_a2*tf(i)^2 + j2_a3*tf(i)^3;
    xx3(i:1) = j1_a1 + 2*j1_a2*tf(i) + 3*j1_a3*tf(i)^2;
    xx4(i:1) = j2_a1 + 2*j2_a2*tf(i) + 3*j2_a3*tf(i)^2;
    time = tf(i);
    x = X(i, :)';
    [~,u] = RRB_ode(time, x, K, P, rho, phi);
    U = [U u];
end
size(U);
figure;
% set(0, 'defaultFigureRenderer', 'painters')
subplot(2,3,1)
plot(tf,(X(:,1)),'r')
hold on
plot(tf,xx1,'--b')
grid on;
title('\theta_1 Simulated Trajectory Tracking')

subplot(2,3,2)
plot(tf,(X(:,3)),'r')
hold on
plot(tf,xx3,'b')
grid on;
title('\theta_1 Simulated Velocity Tracking')

subplot(2,3,3)
plot(tf,U(1,:),'b')
grid on;
title('u_1 Control Input')

subplot(2,3,4)
plot(tf,(X(:,2)),'r')
hold on
plot(tf,xx2,'--b')
grid on;
title('\theta_2 Simulated Trajectory Tracking')

subplot(2,3,5)
plot(tf,(X(:,4)),'r')
hold on
plot(tf,xx4,'b')
grid on;
title('\theta_2 Simulated Velocity Tracking')


subplot(2,3,6)
plot(tf,U(2,:),'b')
grid on;
title('u_2 Control Input')
