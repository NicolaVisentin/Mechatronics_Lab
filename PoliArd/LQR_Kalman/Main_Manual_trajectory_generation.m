clear global
clear
close all
clc


% Data

m=0.264;            % mass of the ball
r=0.02;             % radius of the ball
g=9.81;             % gravity
J=2/5*m*r^2;        % inertia of the ball
b=0.0245;           % arm of the motor
lx=0.182;           % longer edge HALF length
ly=0.150;           % shorter edge HALF length

Cx=m*g*b/lx/(m+J/r^2);
Cy=m*g*b/ly/(m+J/r^2);

% Parametri traiettoria (a velocità costante)

t_sim=30;   % durata traiettoria
Nsec=10;    % risoluzione: punti al secondo

% Acquisisci traiettoria

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

t_acquired=linspace(0,t_sim,length(x_points_acquired));
t_plan=linspace(0,t_sim,Nsec*t_sim)';
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

% Calcola controllo e tutto

Dt=1/Nsec;
xp_plan=gradient(x_plan,Dt);
yp_plan=gradient(y_plan,Dt);
xpp_plan=gradient(xp_plan,Dt);
ypp_plan=gradient(yp_plan,Dt);

ux_plan=asin(xpp_plan/Cx);
uy_plan=asin(ypp_plan/Cy);

figure
subplot(2,1,1)
plot(t_plan,rad2deg(ux_plan))
xlabel('t [s]','Interpreter','latex')
ylabel('$\alpha_x$ [deg]','Interpreter','latex')
grid on
box on
title('planned control')
axis tight
subplot(2,1,2)
plot(t_plan,rad2deg(uy_plan))
xlabel('t [s]','Interpreter','latex')
ylabel('$\alpha_y$ [deg]','Interpreter','latex')
grid on
box on
title('planned control')
axis tight

% Save data

track.x=x_plan;
track.y=y_plan;
track.xp=xp_plan;
track.yp=yp_plan;
track.t=t_plan;
track.ux=ux_plan;
track.uy=uy_plan;

save('track_traj','track')

%% Callback function for mouse click

function mouseClickCallback(~,~)

    % Definisci variabili globali per mantenere i punti tra le chiamate

    global x_points y_points

    % Se le variabili non sono ancora state inizializzate, inizializzale

    if isempty(x_points)
        x_points = [];
        y_points = [];
    end

    % Ottieni il punto dove è stato cliccato

    p = get(gca, 'Currentpoint');  % Ottieni la posizione del click rispetto agli assi
    x = p(1,1);  % La posizione x
    y = p(1,2);  % La posizione y

    % Aggiungi il punto all'elenco

    x_points = [x_points, x];
    y_points = [y_points, y];

    % Traccia i punti

    plot(x_points, y_points, 'ro-', 'MarkerFaceColor','r');

    % Verifica se è stato cliccato con il tasto destro

    if strcmp(get(gcf, 'SelectionType'), 'alt')  % 'alt' è il tasto destro

        % Salva la traiettoria nel workspace
        assignin('base', 'x_points_acquired', x_points);
        assignin('base', 'y_points_acquired', y_points);

        % Disabilita la callback per fermare l'esecuzione
        set(gcf, 'WindowButtonDownFcn', []);
        uiresume;

    end

end