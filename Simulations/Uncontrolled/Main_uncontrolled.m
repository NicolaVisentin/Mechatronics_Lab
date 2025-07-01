clear
close all
clc

%% Data and parameters

% Data of the mechanical system

m=0.264;            % mass of the ball
r=0.02;             % radius of the ball
g=9.81;             % gravity
J=2/5*m*r^2;        % inertia of the ball
b=0.0245;           % arm of the motor
lx=0.182;           % longer edge HALF length
ly=0.150;           % shorter edge HALF length

data_sys=[m r g J b lx ly];

%% Theoretical analysis (LINEARISED SYSTEM)

% Plant transfer functions

s=tf('s');

Cx=m*g*b/(m+J/r^2)/lx;
Cy=m*g*b/(m+J/r^2)/ly;

G_x=Cx/s^2;
G_y=Cy/s^2;

% Poles of the uncontrolled plant

poles_x_un=pole(G_x);
poles_y_un=pole(G_y);

figure

subplot(2,1,1)
plot(real(poles_x_un),imag(poles_x_un),'x','LineWidth',2,'MarkerSize',7)
grid on
box on
xline(0)
yline(0)
xlabel('Re [rad/s]','Interpreter','latex')
ylabel('Im [rad/s]','Interpreter','latex')
title('Unontrolled poles (x direction)')
ax=axis; 
dax=0.2*(ax(2)-ax(1));
day=0.2*(ax(4)-ax(3));
axis([ax(1)-dax ax(2)+dax ax(3)-day ax(4)+day])

subplot(2,1,2)
plot(real(poles_y_un),imag(poles_y_un),'x','LineWidth',2,'MarkerSize',7)
grid on
box on
xline(0)
yline(0)
xlabel('Re [rad/s]','Interpreter','latex')
ylabel('Im [rad/s]','Interpreter','latex')
title('Unontrolled poles (y direction)')
ax=axis; 
dax=0.2*(ax(2)-ax(1));
day=0.2*(ax(4)-ax(3));
axis([ax(1)-dax ax(2)+dax ax(3)-day ax(4)+day])

% Bode diagram of the uncontrolled plant

figure
bode(G_x)
grid on
box on
title('Uncontrolled Bode plot (x direction)')

figure
bode(G_y)
grid on
box on
title('Uncontrolled Bode plot (y direction)')

[~,pm_un_x,~,wc_un_x]=margin(G_x);  % phase margin and cutting frequency
[~,pm_un_y,~,wc_un_y]=margin(G_y);  % phase margin and cutting frequency

fprintf('\nphase margin (x): %.3f deg\ncutting frequency (x): %.3f rad/s \nphase margin (y): %.3f deg \ncutting frquency (y): %3.f rad/s\n\n', pm_un_x, wc_un_x, pm_un_y, wc_un_y)

%% Simulink simulation (NONLINEAR SYSTEM)

% Simulation time

t_sim=20;
    
% Initial conditions

x0=0.05;
xp0=0;

y0=0;
yp0=0;

zx0=[x0 xp0]';
zy0=[y0 yp0]'; 

    % Control
    
    tk=linspace(0,20,1001)';
    theta=0:0.1:2*pi;
    thetap=2*pi/5;
    
    ux= asin((-0.05*cos(theta)*thetap^2)/Cx);
    uy= asin((-0.05*sin(theta)*thetap^2)/Cy);
    
    ux_k=interp1(theta,ux,tk); ux_k(312:end)=[];
    uy_k=interp1(theta,uy,tk); uy_k(312:end)=[];
    tk(312:end)=[];

    load('track_traj.mat')
    tk=track.t;
    ux_k=track.ux;
    uy_k=track.uy;

% simulate

out=sim('plant_uncontrolled.slx');

%% Results

t=out.tout;

% angles of the motors (alpha) and inclination of the platform (theta)

alpha_x=squeeze(out.alpha_x);
alpha_y=squeeze(out.alpha_y);
theta_x=squeeze(out.theta_x);
theta_y=squeeze(out.theta_y);

figure
subplot(2,1,1)
title('Motor angles')
hold on
plot(t,rad2deg(alpha_x))
plot(t,rad2deg(alpha_y))
grid on
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
xlabel('x [mm]','Interpreter','latex')
ylabel('y [mm]','Interpreter','latex')
title('Trajectory of the ball')
axis([-1000*1.1*lx 1000*1.1*lx -1000*1.1*ly 1000*1.1*ly])
hold off

% animation

filename='animation.gif';
FigTag=figure;
for ii=1:100:length(t)

    clf(FigTag)
    
    x_ii=x(ii);
    y_ii=y(ii);

    hold on
    grid on
    plot(1000*x(1:ii),1000*y(1:ii),'k--')
    rectangle('Position',[-1000*lx -1000*ly 1000*2*lx 1000*2*ly],'LineWidth',1.5);
    hp=scatter(1000*x_ii,1000*y_ii,'fill');
    hp.MarkerEdgeColor='k';
    hp.MarkerFaceColor='k';
    hp.LineWidth=4;
    xlabel('$x\;[mm]$','Interpreter','LaTex')
    ylabel('$y\;[mm]$','Interpreter','LaTex')
    title('Animation')
    axis([-1000*1.1*lx 1000*1.1*lx -1000*1.1*ly 1000*1.1*ly])
    box on
    drawnow

    % save animation
    frame=getframe(FigTag);  % [ottieni il frame corrente della figura]
    img=frame2im(frame);     % [converti il frame in immagine RGB]
    [imind,cm]=rgb2ind(img,256);  % [converti l'immagine in formato indicizzato con 256 colori]
    if ii==1
        % [crea il file GIF; imposta quante volte la vuoi loopare e il ritardo tra i loop]
        imwrite(imind, cm, filename, 'gif', 'Loopcount', 1, 'DelayTime', 1);
    else
        % [aggiungi i frame successivi al file GIF e imposta il frame rate]
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/30);
    end

end