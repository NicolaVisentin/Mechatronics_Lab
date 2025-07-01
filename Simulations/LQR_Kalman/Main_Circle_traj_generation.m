clear 
close all
clc

%% Data

m=0.264;            % mass of the ball
r=0.02;             % radius of the ball
g=9.81;             % gravity
J=2/5*m*r^2;        % inertia of the ball
b=0.0245;           % arm of the motor
lx=0.182;           % longer edge HALF length
ly=0.150;           % shorter edge HALF length

Cx=m*g*b/lx/(m+J/r^2);
Cy=m*g*b/ly/(m+J/r^2);

%% Parameters

r=0.05;   % radius
T=5;      % period
N=5;      % number of periods
res=0.1;  % temporal resolution (seconds)

%% Generation

% Start in (r,0) (first 4 seconds)

t=0:res:4;
x=r*ones(1,length(t));
y=zeros(1,length(t));
xp=zeros(1,length(t));
yp=zeros(1,length(t));
ux=zeros(1,length(t));
uy=zeros(1,length(t));
theta=zeros(1,length(t));
thetap=zeros(1,length(t));
thetapp=zeros(1,length(t));

%%%% Triangular acceleration profile (max acc more or less =0.866, final velocity more or less =2*pi/5 --> t=3 sec) %%%%

    % accelerate from thetapp=0, with linear acceleration, in 1.5 sec (more or less 2*pi/(5*0.866))
    
    t_seg=res:res:1.5;
    t=[t t(end)+t_seg];
    
    thetappp_seg=0.866*ones(1,length(t_seg));
    thetapp_seg=thetappp_seg.*t_seg+0;
    thetap_seg=0.5*thetappp_seg.*t_seg.^2+0*t_seg+0;
    theta_seg=1/6*thetappp_seg.*t_seg.^3+0.5*0*t_seg.^2+0*t_seg+0;
    
    x=[x r*cos(theta_seg)];
    y=[y r*sin(theta_seg)];
    
    xp=[xp -r*sin(theta_seg).*thetap_seg];
    yp=[yp r*cos(theta_seg).*thetap_seg];
    
    xpp=-r*cos(theta_seg).*thetap_seg.^2-r*sin(theta_seg).*thetapp_seg;
    ypp=-r*sin(theta_seg).*thetap_seg.^2+r*cos(theta_seg).*thetapp_seg;
    
    ux=[ux asin(xpp/Cx)];
    uy=[uy asin(ypp/Cy)];

    theta=[theta theta_seg];
    thetap=[thetap thetap_seg];
    thetapp=[thetapp thetapp_seg];

    % decelerate until thetapp=0 and speed = 2*pi/5

    t_seg=res:res:1.5;
    t=[t t(end)+t_seg];
    
    thetappp_seg=-0.866*ones(1,length(t_seg));
    thetapp_seg=thetappp_seg.*t_seg+thetapp(end);
    thetap_seg=0.5*thetappp_seg.*t_seg.^2+thetapp(end)*t_seg+thetap(end);
    theta_seg=1/6*thetappp_seg.*t_seg.^3+0.5*thetapp(end)*t_seg.^2+thetap(end)*t_seg+theta(end);
    
    x=[x r*cos(theta_seg)];
    y=[y r*sin(theta_seg)];
    
    xp=[xp -r*sin(theta_seg).*thetap_seg];
    yp=[yp r*cos(theta_seg).*thetap_seg];
    
    xpp=-r*cos(theta_seg).*thetap_seg.^2-r*sin(theta_seg).*thetapp_seg;
    ypp=-r*sin(theta_seg).*thetap_seg.^2+r*cos(theta_seg).*thetapp_seg;
    
    ux=[ux asin(xpp/Cx)];
    uy=[uy asin(ypp/Cy)];

    theta=[theta theta_seg];
    thetap=[thetap thetap_seg];
    thetapp=[thetapp thetapp_seg];

% Start circling around with constant speed

t_seg=res:res:T*N;
t=[t t(end)+t_seg];

thetapp_seg=zeros(1,length(t_seg));
thetap_seg=thetap(end)*ones(1,length(t_seg));
theta_seg=0.5*thetapp_seg.*t_seg.^2+thetap_seg.*t_seg+theta(end);

x=[x r*cos(theta_seg)];
y=[y r*sin(theta_seg)];

xp=[xp -r*sin(theta_seg).*thetap_seg];
yp=[yp r*cos(theta_seg).*thetap_seg];

xpp=-r*cos(theta_seg).*thetap_seg.^2-r*sin(theta_seg).*thetapp_seg;
ypp=-r*sin(theta_seg).*thetap_seg.^2+r*cos(theta_seg).*thetapp_seg;

ux=[ux asin(xpp/Cx)];
uy=[uy asin(ypp/Cy)];

theta=[theta theta_seg];
thetap=[thetap thetap_seg];
thetapp=[thetapp thetapp_seg];

%%%

% figure
% plot(nan,nan)
% axis([-0.1 0.1 -0.1 0.1])
% grid on
% box on
% 
% for ii=1:length(t)
% 
%    clf
%    plot(x(ii),y(ii),'ro')
%    axis([-0.1 0.1 -0.1 0.1])
%    grid on
%    box on
%    drawnow
% 
% end
% 
% figure
% plot(nan,nan)
% axis([0 30 -0.1 0.1])
% grid on
% box on
% 
% for ii=1:length(t)
% 
%    clf
%    plot(t(ii),ux(ii),'ro')
%    axis([0 30 -0.1 0.1])
%    grid on
%    box on
%    drawnow
% 
% end

figure
plot(t,x)
grid on
title('x')

figure
plot(t,y)
grid on
title('y')

figure
plot(t,ux)
grid on
title('ux')

figure
plot(t,uy)
grid on
title('uy')

figure
plot(t,xp)
grid on
title('xp')

figure
plot(t,yp)
grid on
title('yp')

figure
plot(t,theta)
grid on
title('theta')

figure
plot(t,thetap)
grid on
title('thetap')

figure
plot(t,thetapp)
grid on
title('thetapp')

%% Save

circle.x=x';
circle.y=y';
circle.xp=xp';
circle.yp=yp';
circle.t=t';
circle.ux=ux';
circle.uy=uy';

save('circle_traj','circle');
