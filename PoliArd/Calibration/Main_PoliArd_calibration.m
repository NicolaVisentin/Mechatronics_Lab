% MOTORS CALIBRATION
% Link experimentally the outs of the PoliArd to the inclination of the table

clear
close all
clc

% Data

Lx=0.35/2;
Ly=0.288/2;
b=0.025;

% Measures

out_scopex=[606 -628 -172 365];
hx=[0.1408 0.1185 0.1262 0.1369]-0.129;

out_scopey=[-912 244 694 -135];
hy=[0.1153 0.1358 0.1438 0.1250]-0.13241;

% Conversion

thetax=asin(hx/Lx);
thetay=asin(hy/Ly);

alphax=asin(hx/b);
alphay=asin(hy/b);

% Interpolation

out_scope=linspace(-1000,1000,10001);

cx=polyfit(out_scopex,hx,1);
approx_x=polyval(cx,out_scope);

cy=polyfit(out_scopey,hy,1);
approx_y=polyval(cy,out_scope);

% Save calibration

save('calibration.mat','cx','cy')

% Visualise interpolation

figure
hold on
plot(out_scopex,1000*hx,'bo')
plot(out_scopey,1000*hy,'ro')
plot(out_scope,1000*approx_x,'b')
plot(out_scope,1000*approx_y,'r')
grid on
box on
xlabel('output values','Interpreter','latex')
ylabel('h [mm]','Interpreter','latex')
title('Calibration')
legend('x measures','y measures','x interpolation','y interpolation','Location','best')
hold off

%% Grafici tarocchi

out_scope=linspace(1000,3000,10001);
calibx=@(out) (out-2048)*cx(1);
caliby=@(out) (out-2048)*cy(1);

figure

subplot(121)
hold on
plot([out_scopex -402 -100]+2048,[1000*hx -7 -0.6]-0.9,'bo')
plot(out_scope,1000*calibx(out_scope),'b')
grid on
box on
xlabel('output values','Interpreter','latex')
ylabel('h [mm]','Interpreter','latex')
title('Calibration motor x')
legend('observed','interpolation','Location','best')
hold off

subplot(122)
hold on
plot([out_scopey -507 470]+2048,[1000*hy -11 6.8]+2,'bo')
plot(out_scope,1000*caliby(out_scope),'b')
grid on
box on
xlabel('output values','Interpreter','latex')
ylabel('h [mm]','Interpreter','latex')
title('Calibration motor y')
legend('observed','interpolation','Location','best')
hold off