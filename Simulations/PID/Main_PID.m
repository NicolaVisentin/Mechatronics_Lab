clear
%close all
clc

%% Data

m=0.264;            % mass of the ball
r=0.02;             % radius of the ball
g=9.81;             % gravity
J=2/5*m*r^2;        % inertia of the ball 
b=0.0245;           % arm of the motor
lx=0.182;           % longer edge HALF length
ly=0.150;           % shorter edge HALF length

data_sys=[m r g J b lx ly];

% Plant transfer functions

s=tf('s');

Cx=m*g*b/(m+J/r^2)/lx;
Cy=m*g*b/(m+J/r^2)/ly;

G_x=Cx/s^2;
G_y=Cy/s^2;

% State space formulation

Ax=[0 1; 0 0];
Bx=[0 Cx]';
CCx=[1 0];
Dx=0;

Ay=[0 1; 0 0];
By=[0 Cy]';
CCy=[1 0];
Dy=0;

%% PID controller design

% Controller gains

Kp_x=5;
Ki_x=5;
Kd_x=3;

Kp_y=Kp_x/1.2;
Ki_y=Ki_x/1.2;
Kd_y=Kd_x/1.2;

% High frequency pole

N_x=15;   % add high freq pole for real implementation
N_y=15;   % add high freq pole for real implementation

% Input limitations [deg]

ux_lim=45;
uy_lim=45;

% Anti wind-up parameter (if used)

Kb_x=0;
Kb_y=0;

% PID controller transfer function

R_x=Kp_x+Kd_x*s*(N_x/(s+N_x))+Ki_x/s;
R_y=Kp_y+Kd_y*s*(N_y/(s+N_y))+Ki_y/s;

% Open loop and closed loop tranfer functions

L_x=R_x*G_x;
L_y=R_y*G_y;

F_x=feedback(L_x,1);
F_y=feedback(L_y,1);

% Bode plots

figure
bode(L_x)
grid on
box on
title('Open loop Bode plot (x direction)')

figure
bode(L_y)
grid on
box on
title('Open loop Bode plot (y direction)')

[~,pm_x,~,wc_x]=margin(L_x);  % phase margin
[~,pm_y,~,wc_y]=margin(L_y);  % phase margin

fprintf('X DIRECTION: PERFORMANCE\n\tphase margin: %.3f deg\n\tcutting frequency: %.3f rad/s \nY DIRECTION: PERFORMANCE\n\tphase margin: %.3f deg \n\tcutting frquency: %3.f rad/s\n\n', pm_x, wc_x, pm_y, wc_y)

% Step response

fprintf('X DIRECTION: PERFORMANCE\n')
infox=stepinfo(F_x);
disp(infox)
figure
step(F_x)
grid on
title('Step response - x direction')

fprintf('Y DIRECTION: PERFORMANCE\n')
infoy=stepinfo(F_y);
disp(infoy)
figure
step(F_y)
grid on
title('Step response - y direction')

%% Simulink simulation

% Simulation time

t_sim=10;

% Initial conditions

x0=0;
y0=0;

xp0=0;
yp0=0;

zx0=[x0 xp0]';
zy0=[y0 yp0]';

% Reference

x_ref=0.05;
y_ref=0.05;

% Simulation

out=sim('PID.slx');

%% Results

t=out.tout;

% angles of the motors (alpha) and inclination of the platform (theta)

alpha_x=out.alpha_x;
alpha_y=out.alpha_y;
theta_x=out.theta_x;
theta_y=out.theta_y;

figure
subplot(2,1,1)
title('Motor angles')
hold on
plot(t,rad2deg(alpha_x))
plot(t,rad2deg(alpha_y))
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\alpha$ [deg]','Interpreter','latex')
legend('\alpha_x','\alpha_y','Location','best')
axis tight
hold off
subplot(2,1,2)
title('Platform angles')
hold on
plot(t,rad2deg(theta_x))
plot(t,rad2deg(theta_y))
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\theta$ [deg]','Interpreter','latex')
legend('\theta_x','\theta_y','Location','best')
axis tight
hold off

% position of the ball

x=squeeze(squeeze(out.x));
y=squeeze(squeeze(out.y));

figure
hold on
plot(t,1000*x)
plot(t,1000*y)
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('position [mm]','Interpreter','latex')
title('Position of the ball')
legend('x','y','Location','best')
axis tight
hold off

figure
hold on
plot(1000*x,1000*y,'b')
plot(1000*zx0(1), 1000*zy0(2),'bo')
rectangle('Position',[-1000*lx -1000*ly 1000*2*lx 1000*2*ly],'LineWidth',1.5)
grid on
box on
xlabel('x [mm]','Interpreter','latex')
ylabel('y [mm]','Interpreter','latex')
title('Trajectory of the ball')
axis([-1000*1.1*lx 1000*1.1*lx -1000*1.1*ly 1000*1.1*ly])
hold off

% % animation
% 
% filename='animation.gif';
% FigTag=figure;
% for ii=1:100:length(t)
% 
%     clf(FigTag)
% 
%     x_ii=x(ii);
%     y_ii=y(ii);
% 
%     hold on
%     grid on
%     plot(1000*x(1:ii),1000*y(1:ii),'k--')
%     rectangle('Position',[-1000*lx -1000*ly 1000*2*lx 1000*2*ly],'LineWidth',1.5);
%     hp=scatter(1000*x_ii,1000*y_ii,'fill');
%     hp.MarkerEdgeColor='k';
%     hp.MarkerFaceColor='k';
%     hp.LineWidth=4;
%     xlabel('$x\;[mm]$','Interpreter','LaTex')
%     ylabel('$y\;[mm]$','Interpreter','LaTex')
%     title('Animation')
%     axis([-1000*1.1*lx 1000*1.1*lx -1000*1.1*ly 1000*1.1*ly])
%     box on
%     drawnow
% 
%     % save animation
%     frame=getframe(FigTag);  % [ottieni il frame corrente della figura]
%     img=frame2im(frame);     % [converti il frame in immagine RGB]
%     [imind,cm]=rgb2ind(img,256);  % [converti l'immagine in formato indicizzato con 256 colori]
%     if ii==1
%         % [crea il file GIF; imposta quante volte la vuoi loopare e il ritardo tra i loop]
%         imwrite(imind, cm, filename, 'gif', 'Loopcount', 1, 'DelayTime', 1);
%     else
%         % [aggiungi i frame successivi al file GIF e imposta il frame rate]
%         imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/30);
%     end
% 
% end