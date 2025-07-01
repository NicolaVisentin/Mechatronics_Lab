clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   NOTATION                                                              %
%       - zx,zy:        [x xp]', [y yp]' (states)                         %
%       - u_x,u_y:      alphax, alphay (control)                          %
%       - P:            final state weighting matrix                      %
%       - PP:           Riccati matrix                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data and parameters definition

m=0.264;            % mass of the ball
r=0.02;             % radius of the ball
g=9.81;             % gravity
J=2/5*m*r^2;        % inertia of the ball
b=0.025;           % arm of the motor
lx=0.35/2;           % longer edge HALF length
ly=0.288/2;           % shorter edge HALF length

Cx=m*g*b/lx/(m+J/r^2);
Cy=m*g*b/ly/(m+J/r^2);

data_sys=[m r g J b lx ly];

%% Infinite time LQR

% weights on delta_x and delta_u; final target; integrative term

r_delta_alphax = 0.05/(1)^2;
r_delta_alphay = 0.05/(1)^2;
q_delta_x = 1/(1)^2;
q_delta_y = 1/(1)^2;
q_delta_xp = 0;
q_delta_yp = 0;
q_wx=.000001;
q_wy=.000001;

Rx_lqr=r_delta_alphax;                  % L=0.5*(u'*R*u+x'*Q*x)
Qx_lqr=[q_delta_x 0; 0 q_delta_xp];     % L=0.5*(u'*R*u+x'*Q*x)

Ry_lqr=r_delta_alphay;                  % L=0.5*(u'*R*u+x'*Q*x)
Qy_lqr=[q_delta_y 0; 0 q_delta_yp];     % L=0.5*(u'*R*u+x'*Q*x)

zx_f=[0.05 0]';
zy_f=[0.05 0]';
alphax_f=0;
alphay_f=0;

% A, B and C matrices

A=[0 1; 0 0];
Bx=[0; Cx];
By=[0; Cy];
C=[1 0];

poles_un=eig(A);

% Extended matrices (for steady-state error)

A_ex=[A, zeros(2,1); C, 0];
Bx_ex=[Bx; 0];
By_ex=[By; 0];

Qx_ex=[Qx_lqr, zeros(2,1); zeros(1,2), q_wx];
Qy_ex=[Qy_lqr, zeros(2,1); zeros(1,2), q_wy];

% design of the LQR using matlab function

[Kx_lqr,~,poles_conx]=lqr(A_ex,Bx_ex,Qx_ex,Rx_lqr);
[Ky_lqr,~,poles_cony]=lqr(A_ex,By_ex,Qy_ex,Ry_lqr);

% plot poles of the controlled and uncontrolled system

figure
hold on
plot(real(poles_un),imag(poles_un),'bx','LineWidth',2)
plot(real(poles_conx),imag(poles_conx),'rx','LineWidth',2)
xline(0)
yline(0)
grid on
box on
xlabel('Re','Interpreter','latex')
ylabel('Im','Interpreter','latex')
title('Poles (along x or y direction)')
legend('uncontrolled system','controlled system')
hold off

%% Kalman filter design

% Covariance matrices

Qxk=diag([100, 1000]);         % disturbance covariance
Qyk=diag([100, 1000]);         % disturbance covariance
Rxk=.1;          % noise covariance
Ryk=.1;          % noise covariance

% observer gain matrix

[~,Kx_obs,polesx_kalm]=icare(A',C',Qxk,Rxk);
Kx_obs=Kx_obs';

[~,Ky_obs,polesy_kalm]=icare(A',C',Qyk,Ryk);
Ky_obs=Ky_obs';

% observer state matrices

A_obs=A;
Bx_obs=[Bx Kx_obs];
By_obs=[By Ky_obs];
C_obs=eye(2);
D_obs=zeros(2);

% plot all poles

figure
hold on
plot(real(poles_un),imag(poles_un),'bx','LineWidth',2)
plot(real(poles_conx),imag(poles_conx),'rx','LineWidth',2)
plot(real(polesx_kalm),imag(polesx_kalm),'kx','LineWidth',2)
xline(0)
yline(0)
grid on
box on
xlabel('Re [rad/s]','Interpreter','LaTex')
ylabel('Im [rad/s]','Interpreter','LaTex')
legend('Uncontrolled system','Controlled system','Observer','Interpreter','LaTex','Location','Best')
title('Poles')
ax=axis;
axis([1.2*ax(1)+1 1.2*ax(2)+1 1.2*ax(3)+1 1.2*ax(4)+1])
hold off

%% Simulink

% load nominal trajectory

load('circle_traj')
x_plan=circle.x;
y_plan=circle.y;
t_plan=circle.t;

% parameters

zx0=[0.05,0]';          % initial conditions
zy0=[0,0]';          % initial conditions

ICx_obs=zx0-zx_f;    % observer initial conditions (on delta_zx)
ICy_obs=zy0-zy_f;    % observer initial conditions (on delta_zy)

tf=60;

% simulate and extract results

out=sim('LQR_SteadyState_Kalman.slx');

t_simul=squeeze(out.tout);
x_simul=squeeze(out.x_sim);
alphax_simul=squeeze(out.alphax_sim);
x_ref=squeeze(out.x_ref);
alphax_ref=squeeze(out.alphax_ref);
y_simul=squeeze(out.y_sim);
alphay_simul=squeeze(out.alphay_sim);
y_ref=squeeze(out.y_ref);
alphay_ref=squeeze(out.alphay_ref);
delta_zx_hat=out.delta_zx_hat;
delta_zy_hat=out.delta_zy_hat;
xp_simul=squeeze(out.xp_sim);
yp_simul=squeeze(out.yp_sim);

%% Plots

% x, alpha_x and y, alpha_y in time

figure

subplot(2,1,1)
hold on
plot(t_simul,1000*x_simul,'b')
plot(t_simul,1000*x_ref,'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('x [mm]','Interpreter','latex')
title('x position')
legend('$x$','$x_{ref}$','Interpreter','latex','location','best')
hold off

subplot(2,1,2)
hold on
plot(t_simul,rad2deg(alphax_simul),'r')
plot(t_simul,rad2deg(ones(1,length(t_simul))*alphax_ref),'r--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\alpha$ [deg]','Interpreter','latex')
title('Control in x')
legend('$\alpha_x$','$\alpha_{x,ref}$','Interpreter','latex','location','best')
hold off

figure

subplot(2,1,1)
hold on
plot(t_simul,1000*y_simul,'b')
plot(t_simul,1000*y_ref,'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('y [mm]','Interpreter','latex')
title('y position')
legend('$y$','$y_{ref}$','Interpreter','latex','location','best')
hold off

subplot(2,1,2)
hold on
plot(t_simul,rad2deg(alphay_simul),'r')
plot(t_simul,rad2deg(ones(1,length(t_simul))*alphay_ref),'r--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\alpha$ [deg]','Interpreter','latex')
title('Control in y')
legend('$\alpha_y$','$\alpha_{y,ref}$','Interpreter','latex','location','best')
hold off

% trajectory

figure
hold on
plot(1000*x_simul,1000*y_simul,'b')
plot(1000*x_simul(1), 1000*y_simul(1),'bo')
rectangle('Position',[-1000*lx -1000*ly 1000*2*lx 1000*2*ly],'LineWidth',1.5)
grid on
box on
xlabel('x [mm]','Interpreter','latex')
ylabel('y [mm]','Interpreter','latex')
title('Trajectory of the ball')
axis([-1000*1.1*lx 1000*1.1*lx -1000*1.1*ly 1000*1.1*ly])
hold off

% observer

figure

subplot(2,1,1)
hold on
plot(t_simul,1000*x_simul,'b')
plot(t_simul,1000*x_ref,'b--')
plot(t_simul,1000*(delta_zx_hat(:,1)+x_ref),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('x [mm]','Interpreter','latex')
title('x position')
legend('$x$','$x_{ref}$','$\hat{x}$','Interpreter','latex','location','best')
hold off

subplot(2,1,2)
hold on
plot(t_simul,xp_simul,'b')
plot(t_simul,ones(length(t_simul),1)*zx_f(2),'b--')
plot(t_simul,(delta_zx_hat(:,2)+zx_f(2)),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\dot{x}$ [m/s]','Interpreter','latex')
title('x velocity')
legend('$\dot{x}$','$\dot{x}_{ref}$','$\dot{\hat{x}}$','Interpreter','latex','location','best')
hold off

figure

subplot(2,1,1)
hold on
plot(t_simul,1000*y_simul,'b')
plot(t_simul,1000*y_ref,'b--')
plot(t_simul,1000*(delta_zy_hat(:,1)+y_ref),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('y [mm]','Interpreter','latex')
title('y position')
legend('$y$','$y_{ref}$','$\hat{y}$','Interpreter','latex','location','best')
hold off

subplot(2,1,2)
hold on
plot(t_simul,yp_simul,'b')
plot(t_simul,ones(length(t_simul),1)*zy_f(2),'b--')
plot(t_simul,(delta_zy_hat(:,2)+zy_f(2)),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\dot{y}$ [m/s]','Interpreter','latex')
title('y velocity')
legend('$\dot{y}$','$\dot{y}_{ref}$','$\dot{\hat{y}}$','Interpreter','latex','location','best')
hold off

% animation

filename='animation.gif';
FigTag=figure;
for ii=1:100:length(t_simul)

    clf(FigTag)
    
    x_ii=x_simul(ii);
    y_ii=y_simul(ii);

    hold on
    grid on
    plot(1000*x_simul(1:ii),1000*y_simul(1:ii),'k--')
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