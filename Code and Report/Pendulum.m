% =========================================================================
%                                FINAL WORK                               
% Course: Modeling and Quantification of Uncertainties in Physical Systems
% Author: Vítor Marchiori
% Date: June 2025
% System: Simple Pendulum with Uncertain Length
% Objective: Model and quantify the uncertainty in the pendulum length (L)
% =========================================================================

clear; clc; close all;

% ----------------------------------------------------------------------- %
%                           TIME DISCRETIZATION                           %
% ----------------------------------------------------------------------- %

t0 = 0.0;                       % Start Time [s]
t1 = 30.0;                      % End Time [s]
Ndt = 300;                      % Number of Time Steps
tspan = linspace(t0, t1, Ndt);  % Time Vector

% ----------------------------------------------------------------------- %
%                               PARAMETERS                                %
% ----------------------------------------------------------------------- %

g = 9.81;            % Gravity [m/s²]
c = 0.1;             % Damping Coefficient [N·s/m]
m = 1.0;             % Pendulum Mass [kg]

theta0 = 0.2;        % Initial Angular Displacement [rad]
omega0 = 0.0;        % Initial Angular Velocity [rad/s]

% ----------------------------------------------------------------------- %
%                          MONTE CARLO SETTINGS                           %
% ----------------------------------------------------------------------- %

Ns = 1024;           % Number of Samples
mean_L = 1.0;        % Mean Length [m]
coefvar_L = 0.1;     % Coeficiente de variação de L (10%)

% Generate Length Samples from Gamma Distribution
rng(30081984);             % Seed for Reproducibility
L = gamrnd(1/coefvar_L^2, mean_L*coefvar_L^2, [Ns,1]);

% Allocate memory
Qtheta = zeros(Ndt, Ns);    % Angle
Qomega = zeros(Ndt, Ns);    % Angular Velocity

% ----------------------------------------------------------------------- %
%                       MONTE CARLO SIMULATION LOOP                       %
% ----------------------------------------------------------------------- %

for n = 1:Ns
    IC = [theta0; omega0];

    % Sample Dependent Parameter of L
    w_n = sqrt(g / L(n));
    zeta = c / (2 * m * L(n)^2 * w_n);

    % Differential Equations of the Damped Linearized Pendulum
    dydt = @(t,y)[y(2); -2*zeta*w_n*y(2) - w_n^2*y(1)];

    [~, Y] = ode45(dydt, tspan, IC);

    Qtheta(:,n) = Y(:,1);  % Angle
    Qomega(:,n) = Y(:,2);  % Angular Velocity
end

% ----------------------------------------------------------------------- %
%                       STATISTICAL POST-PROCESSING                       %
% ----------------------------------------------------------------------- %

theta_mean = mean(Qtheta, 2);
theta_std  = std(Qtheta, 0, 2);
theta_low  = prctile(Qtheta, 2.5, 2);
theta_upp  = prctile(Qtheta, 97.5, 2);
Qtheta_tmean = mean(Qtheta);

% ----------------------------------------------------------------------- %
%                                  PLOTS                                  %
% ----------------------------------------------------------------------- %

figure('Color', 'w', 'Position', [100, 100, 800, 500]);
fill([tspan fliplr(tspan)], [theta_upp' fliplr(theta_low')], [1 1 0]*0.9, ...
    'EdgeColor','none','FaceAlpha',0.3); hold on;
plot(tspan, theta_mean, 'b-', 'LineWidth', 2);
plot(tspan, theta_mean + theta_std, 'r--', 'LineWidth', 2);
plot(tspan, theta_mean - theta_std, 'r--', 'LineWidth', 2);
legend('Intervalo de Confian\c{c}a','M\''{e}dia','M\''{e}dia $\pm$ Desvio Padr{\~a}o', 'Interpreter', 'latex');
xlabel('Tempo [s]', 'FontSize', 16, 'Interpreter', 'latex'); 
ylabel('$\theta$ [rad]', 'FontSize', 16, 'Interpreter', 'latex'); 
%title('Resposta Média do Pêndulo com Banda de Confiança 95%', 'FontSize', 18, 'FontWeight', 'bold'); 
grid on; box on; 
ax = gca;
ax.FontSize = 16; 
ax.XMinorGrid = 'on'; 
ax.YMinorGrid = 'on'; 
ax.GridAlpha = 0.3; 

figure('Color', 'w', 'Position', [100, 100, 800, 500]);
histogram(Qtheta_tmean, 'Normalization','pdf','FaceAlpha', 0.5,  ...
    'FaceColor', [0 0 1], 'EdgeColor', [0 0 0]);
hold on;
[f, xi] = ksdensity(Qtheta_tmean);
plot(xi, f, 'r-', 'LineWidth', 2);
%legend('Histograma','KDE');
xlabel('$\theta$', 'FontSize', 16, 'Interpreter', 'latex'); 
ylabel('Densidade de Probabilidade', 'FontSize', 16, 'Interpreter', 'latex'); 
%title('Histograma e KDE da Média Temporal de $\theta$', 'FontSize', 18, 'FontWeight', 'bold'); 
grid on; box on; 
ax = gca;
ax.FontSize = 16; 
ax.XMinorGrid = 'on'; 
ax.YMinorGrid = 'on'; 
ax.GridAlpha = 0.3; 

% ----------------------------------------------------------------------- %
%                               ANIMATION                                 %
% ----------------------------------------------------------------------- %

%output_folder = 'frames_pendulo';
%if ~exist(output_folder, 'dir')
%    mkdir(output_folder);
%end

idx = randi(Ns);                
theta_anim = Qtheta(:, idx);    
L_anim = L(idx);                

figure('Color','w', 'Position', [100, 100, 800, 500]);
for k = 1:length(tspan)

    x =  L_anim * sin(theta_anim(k));
    y = -L_anim * cos(theta_anim(k));

    clf;
    plot([0 x], [0 y], 'k-', 'LineWidth', 2); hold on;
    plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); 

    axis equal;
    xlim([-1.2 1.2]*mean_L);
    ylim([-1.2 0.2]*mean_L);
    %title(sprintf('Amostra #%d - t = %.2f s', idx, tspan(k)), 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('$x$ [m]', 'FontSize', 16, 'Interpreter', 'latex'); 
    ylabel('$y$ [m]', 'FontSize', 16, 'Interpreter', 'latex'); 
    grid on; box on; 
    ax = gca;
    ax.FontSize = 16; 
    ax.XMinorGrid = 'on'; 
    ax.YMinorGrid = 'on'; 
    ax.GridAlpha = 0.3; 
    drawnow;

    %filename = sprintf('%s/frame_%03d.png', output_folder, k);
    %exportgraphics(gcf, filename, 'BackgroundColor', 'white');

end