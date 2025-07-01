clear global
clear
close all
clc

%% Data

b=0.0245;           % arm of the motor
lx=0.182;           % longer edge HALF length
ly=0.150;           % shorter edge HALF length

load('calibration.mat')

%% Simulink

% Circle trajectory (internal, simulink generated)

T_circle=2;          % period of the circle

% Manual-defined nominal trajectory

track_duration=50;   % if user defined trajectory

Fig=figure;
axis([-lx lx -ly ly])
hold on
grid on
box on
xlabel('x [m]','Interpreter','latex')
ylabel('y [m]','Interpreter','latex')
title('Track your trajectory (right click to stop)')

set(gcf, 'WindowButtonDownFcn', @mouseClickCallback);
uiwait;

t_acquired=linspace(0,track_duration,length(x_points_acquired));
    %t_acquired=linspace(0,5,length(x_points_acquired));
t_plan=linspace(0,track_duration,1*track_duration)';
    %t_plan=0:0.1:5;
x_plan=spline(t_acquired,x_points_acquired,t_plan);
y_plan=spline(t_acquired,y_points_acquired,t_plan);

figure
hold on
plot(x_points_acquired,y_points_acquired,'ro')
plot(x_plan,y_plan,'b')
grid on
box on
axis([-lx lx -ly ly])
xlabel('x [m]','Interpreter','latex')
ylabel('y [m]','Interpreter','latex')
title('Planned trajectory')
legend('selected points','generated trajectory','Interpreter','latex','Location','best')
hold off

% Load nominal trajectory

load('circle_traj.mat')
t_plan=circle.t;
x_plan=circle.x;
y_plan=circle.y;