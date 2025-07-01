clear
close all
clc

% Specifica il percorso del file

filePath='2024-12-19_18-26-32.csv';

% Leggi il file CSV

data=readtable(filePath);

% Accedi ai dati per analisi

time_hard=data.hardware_time-data.hardware_time(1);
time_py=data.python_time;
x=data.x;
y=data.y;
x_ref=data.x_ref;
y_ref=data.y_ref;
ux=data.ux;
uy=data.uy;

% Plots

figure

subplot(4,1,1)
hold on
plot(time_py,x,'b')
plot(time_py,x_ref,'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
xlabel('x [mm]','Interpreter','latex')
title('x position')
legend('x','$x_{ref}$','interpreter','latex','Location','best')
axis tight
hold off

subplot(4,1,2)
hold on
plot(time_py,y,'b')
plot(time_py,y_ref,'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
xlabel('y [mm]','Interpreter','latex')
title('y position')
legend('y','$y_{ref}$','interpreter','latex','Location','best')
axis tight
hold off

subplot(4,1,3)
plot(time_hard,ux)
grid on
box on
xlabel('t [s]','Interpreter','latex')
xlabel('$\alpha_x$ [deg]','Interpreter','latex')
title('x control')
axis tight

subplot(4,1,4)
plot(time_hard,uy)
grid on
box on
xlabel('t [s]','Interpreter','latex')
xlabel('$\alpha_y$ [deg]','Interpreter','latex')
title('y control')
axis tight

figure
plot(x,y)