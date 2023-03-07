
j1_a0 = pi;
j1_a1 = 0;
j1_a2 = -(3*pi)/100;
j1_a3 = pi/500;

j2_a0 = pi/2;
j2_a1 = 0;
j2_a2 = -(3*pi)/200;
j2_a3 = pi/1000;

tf = 10;
traj = [deg2rad(200), deg2rad(125), 0, 0];
[t,traj] = ode45(@RRB_ode,[0,tf], traj);

for i  = 1 : length(t)
    traj1(i,1) = j1_a0 + j1_a1 * t(i) + j1_a2 * t(i)^2 + j1_a3 * t(i)^3;
    traj2(i,1) = j2_a0 + j2_a1 * t(i) + j2_a2 * t(i)^2 + j2_a3 * t(i)^3;
    traj3(i,1) = j1_a1 + 2 * j1_a2 * t(i) + 3 * j1_a3 * t(i)^2;
    traj4(i,1) = j2_a1 + 2 * j2_a2 * t(i) + 3 * j2_a3 * t(i)^2;
end

set(0, 'defaultFigureRenderer', 'painters')
subplot(2,2,1)
plot(t,traj1,t,traj(:,1))
title('Time vs Theta 1')
xlabel('Time');
ylabel('theta 1');

subplot(2,2,3)
plot(t,traj3,t,traj(:,3))
title('Velocity(theta1d) Trajectory')
xlabel('Time');
ylabel('Theta1dot Trajectory');

subplot(2,2,2)
plot(t,traj2,t,traj(:,2))
title('Time vs Theta 2')
xlabel('Time');
ylabel('theta 2');

subplot(2,2,4)
plot(t,traj4,t,traj(:,4))
title('Velocity(theta2d) Trajectory')
xlabel('Time');
ylabel('Theta2dot Trajectory');
