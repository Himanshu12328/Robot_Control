rosshutdown;
clear; close; clc;

% ROS Setup
rosinit;
m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;
r1 = 0.45;
r2 = 0.45;
g = 9.81;
I1 = 0.084;
I2 = 0.084;

y = [];
T = [];

A = [0, 0, 1, 0;
     0, 0, 0, 1;
     0, 0, 0, 0;
     0, 0, 0, 0];
B = [0, 0;
     0, 0;
     1, 0;
     0, 1];
lambda = [-0.6 -0.05 -7 -50];

K = place(A, B, lambda);

toq = [];

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);
tic;
i = 1;
plot1 = [];
t = 0;
while(t < 10)
    t = toc;
    % read the joint states
    jointData = receive(JointStates);
    % inspect the "jointData" variable in MATLAB to get familiar with its structure
    % implement your state feedback controller below
%     K1 = [26.2200 -7.6202 6.7730 -4.8545];
%     K2 = [6.8728 0.9502 2.2068 -1.2374];
    theta1 = jointData.Position(1);
    theta2 = jointData.Position(2);
    theta1dot = jointData.Velocity(1);
    theta2dot = jointData.Velocity(2);

    J1_a0 = pi;
    J1_a1 = 0;
    J1_a2 = -(3*pi)/100;
    J1_a3 = pi/500;
    J2_a0 = pi/2;
    J2_a1 = 0;
    J2_a2 = -(3*pi)/200;
    J2_a3 = pi/1000;

    q = [theta1;
         theta2;
         theta1dot;
         theta2dot];
    y = [y q];
    
    q_dot = [J1_a0 + J1_a1 * t + (J1_a2 * (t)^2) + (J1_a3 * (t)^3);
             J1_a1 + 2 * J1_a2 * t + 3 * (J1_a3 * (t)^2);
             J2_a0 + J2_a1 * t + (J2_a2 * (t)^2) + (J2_a3 * (t)^3);
             J2_a1 + 2 * J2_a2 * t + 3 * (J2_a3 * (t)^2)];

    q_ddot = [2 * J1_a2 + 6 * J1_a3 * t; 
              2 * J2_a2 + 6 * J2_a3 * t];

    G = [- g*m2*r2*sin(theta1 + theta2) - g*l1*m1*sin(theta1) - g*m1*r1*sin(theta1);
         -g*m2*r2*sin(theta1 + theta2)];

    C = [-l1*m2*r2*theta2dot*sin(theta2)*(2*theta1dot + theta2dot);
         l1*m2*r2*theta1dot^2*sin(theta2)];
    M = [m1*r1^2 + I1 + I2 + (m2*(2*l1^2 + 4*cos(theta2)*l1*r2 + 2*r2^2))/2, I2 + (m2*(2*r2^2 + 2*l1*cos(theta2)*r2))/2;
         I2 + (m2*(2*r2^2 + 2*l1*cos(theta2)*r2))/2,                               m2*r2^2 + I2];
    
    Virtual_U = M*(-K * (q - q_dot) + q_ddot ) + C + G;
% 
% 
%     tau1.Data = -K1 * [jointData.Position(1,1), jointData.Position(2,1), jointData.Velocity(1,1), jointData.Velocity(2,1)]';
%     tau2.Data = -K2 * [jointData.Position(1,1), jointData.Position(2,1), jointData.Velocity(1,1), jointData.Velocity(2,1)]';
    tau1.Data = Virtual_U(1);
    tau2.Data = Virtual_U(2);
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    % sample the time, joint state values, and calculated torques here to be plotted at the end
    plot1(i,:) = [t jointData.Position(1,1) jointData.Position(2,1) jointData.Velocity(1,1) jointData.Velocity(2,1) tau1.Data tau2.Data];
    i = i + 1;
end

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

figure;
set(0, 'defaultFigureRenderer', 'painters')
subplot(2,2,1);
plot(plot1(:,1),plot1(:,2));
xlabel('Time');
ylabel('theta1');
subplot(2,2,2);
plot(plot1(:,1),plot1(:,3));
xlabel('Time');
ylabel('theta2');
subplot(2,2,3);
plot(plot1(:,1),plot1(:,4));
xlabel('Time');
ylabel('theta1 dot');
subplot(2,2,4);
plot(plot1(:,1),plot1(:,5));
xlabel('Time');
ylabel('theta2 dot');
figure;
subplot(2,1,1);
plot(plot1(:,1),plot1(:,6));
xlabel('Time');
ylabel('tau1');
subplot(2,1,2);
plot(plot1(:,1),plot1(:,7));
xlabel('Time');
ylabel('tau2');


% disconnect from roscore
rosshutdown;
% plot the trajectories

