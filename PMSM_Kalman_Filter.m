%% Engenharia de Controle e Automação
% Instituto Federal Fluminense 2022.1
% Processamento de Sinais - Signals and Systems
% Aluno: Kaique Guimarães Cerqueira
% Prof.: Alexandre Leite
%
% ========================Projeto de Aplicação============================
%   Filtro de Kalman aplicado a controle de Corrente do motor PMSM
% 
% Vectors:
% x = [x1; x2; x3; x4] = [id; iq; w; theta]
% u = [u1; u2] = [Vd; Vq]
% y = [y1; y2] = [Id; Iq]
% 
% 
clear, clc, close all

%% Making the linearized model
% Motor parameters
R = 1.4; % resistência da armadura
Lq = 5.8e-3;
Ld = 6.6e-3;
J = 1.76e-3; % momento de inércia
b = 0.00038818; % Coef. de atrito
lambda_m = 0.1546; % densidade de fluxo magnético
N = 6; % número de espiras

% Simulation Parameters
zero_arr = zeros(3,1); % Array used for state kalman init, mean noises
Ts = 0.00001; % Sampling period
w = 1; % Disturbance standard deviation (omega)
eta = 0.1*1; % Noise standard deviation

%% Operating point problem
theta_e = 2*pi/180; % 2° (em rad)
w_FOC = 100;
T_FOC = 1;
Iqc_range = 0.03;

% Velocity at the end of the startup
x3_op = w_FOC;
% Initializing acceleration and Iq
x3_dot_op_min = 1000;
Iqc_min = -1;

% Local search
% Dúvida: no código diz-se "searching positive values inside possibilities"
% O que significa isto?
% R: É um problema de otimização, utilizado para encontrar o ponto de
% operação. O ponto é definido em id, iq = cte (id_dot, iq_dot = 0), w =
% cte (aceleração=0)
for Iqc = linspace(0, Iqc_range, 2e6)

    x1_op = Iqc*sin(theta_e);
    x2_op = Iqc*cos(theta_e);
%     The real deal:
%     Função custo utilizada para a otimização do mínimo valor de x3_dot_op
%     (w_dot tender a 0)
    x3_dot_op = (3*N/(2*J)) * (x1_op*x2_op*(Ld-Lq) + lambda_m*x2_op) - (x3_op)*b/J;
    
    if (abs(x3_dot_op) < x3_dot_op_min)
        x3_dot_op_min = abs(x3_dot_op);
        Iqc_min = Iqc;
    end
end
fprintf('[Iqc, x3_op] = ')
[Iqc_min, x3_dot_op_min]
% Iqc = Iqcost
Iqc = Iqc_min;

% Operation point state vector:
x1_op = Iqc*sin(theta_e);
x2_op = Iqc*cos(theta_e);
x3_op = w_FOC;
x4_op = theta_e + w_FOC*T_FOC/2;

% Input vector:
u1_op = R*x1_op - N*Lq*x2_op*x3_op;
u2_op = N*Ld*x1_op*x3_op + R*x2_op;
U = [u1_op;u2_op];

%% Linear state-space matrices (4 states):
% x = [x1; x2; x3; x4] = [id; iq; w; theta]
% 
% A32 = 3*N*((Ld-Lq)*x1_op+lambda_m)/(2*J);
% A = [-R/Ld,         N*Lq*x3_op/Ld,      N*Lq*x2_op/Ld,      0;
%      -N*Ld*x3_op/Lq,        -R/Lq,     -Ld*N*x1_op/Lq,      0;
%     3*N*(Ld-Lq)*x2_op/(2*J),  A32,               -b/J,      0;
%           0                     0,                  1,      0];
% clear A32
%       
% B = [1/Ld, 0;
%      0, 1/Lq;
%      0,    0;
%      0,    0];
%  
% C = [1, 0, 0, 0;
%      0, 1, 0, 0];
%   
% D = [0, 0; 0, 0]; 
% 
% % Dinamic noise variance state vector
% Q_vec = ones(4,1)*(w^2); % w² = variance
% % Dinamic noise covariance Kalman matrix
% Q = diag(Q_vec);
% % usar o comando diag para usar o Q e R no filtro de kalman
% % Perturbation impact matrix
% G = eye(4);
% % Measurement noise variance state vector
% R_vec = ones(2,1)*(eta^2);
% % Measurement noise covariance Kalman matrix
% R = diag(R_vec);

%% Simplification - Linear state-space matrices (3 states):
% x = [x1; x2; x3] = [id; iq; w]

A32 = 3*N*((Ld-Lq)*x1_op+lambda_m)/(2*J);
A = [-R/Ld,         N*Lq*x3_op/Ld,      N*Lq*x2_op/Ld;
     -N*Ld*x3_op/Lq,        -R/Lq,     -Ld*N*x1_op/Lq;
    3*N*(Ld-Lq)*x2_op/(2*J),  A32,               -b/J];    
clear A32
      
B = [1/Ld, 0;
     0, 1/Lq;
     0,    0];

% Simplificação: considerar que há a leitura de todos os estados
% C = [1, 0, 0;
%      0, 1, 0];
C = eye(3);
 
% Dinamic noise variance state vector
Q_vec = ones(3,1)*(w^2); % w² = variance
% Dinamic noise covariance Kalman matrix
Q = diag(Q_vec);
% usar o comando diag para usar o Q e R no filtro de kalman
% Perturbation impact matrix
G = eye(3);
% Measurement noise variance state vector
R_vec = ones(3,1)*(eta^2);
% Measurement noise covariance Kalman matrix
R = diag(R_vec);

%% Running Simulation
out = sim('Kalman_sim_R2021a')

%% Analysing Results - Data scraping
kalman.id = out.kalman(:,1)';
kalman.iq = out.kalman(:,2)';
kalman.omega = out.kalman(:,3)';

ideal.id = out.ideal(1,:);
ideal.iq = out.ideal(2,:);
ideal.omega = out.ideal(3,:);

noisy_measure.id = out.noisy(1,:);
noisy_measure.iq = out.noisy(2,:);
noisy_measure.omega = out.noisy(3,:);

noisy_states.id = out.noisy_states(1,:);
noisy_states.iq = out.noisy_states(2,:);
noisy_states.omega = out.noisy_states(3,:);

cov = out.covariance.^0.5;
dp.id = cov(1,:);
dp.iq = cov(2,:);
dp.omega = cov(3,:);

time = out.tout';

clear cov out

%% 2d Plot
% Id and Iq
figure()
    subplot(2,1,1)
    plot(time, kalman.id + dp.id, '--', 'Color','#EDB120')
    hold on
    plot(time, kalman.id - dp.id, '--', 'Color','#EDB120')
    plot(time, kalman.id, 'Color','#D95319')
    plot(time, ideal.id, 'Color','#77AC30')
    plot(time, noisy_measure.id, 'Color', '#0072BD')
    title("Id")
    axis tight

    subplot(2,1,2)
    plot(time, kalman.iq + dp.iq, '--', 'Color','#EDB120')
    hold on
    plot(time, kalman.iq - dp.iq, '--', 'Color','#EDB120')
    plot(time, kalman.iq, 'Color','#D95319')
    plot(time, ideal.iq, 'Color','#77AC30')
    plot(time, noisy_measure.iq, 'Color', '#0072BD')
    title("Iq")
    xlabel('Time (s)')
    axis tight

% Omega
figure()
    plot(time, kalman.omega + dp.omega, '--', 'Color','#EDB120')
    hold on
    plot(time, kalman.omega - dp.omega, '--', 'Color','#EDB120')
    plot(time, kalman.omega, 'Color','#D95319')
    plot(time, ideal.omega, 'Color','#77AC30')
    plot(time, noisy_measure.omega, 'Color', '#0072BD')
    title("Omega(\omega)")
    xlabel('Time (s)')
    ylabel('Probability Density')
    axis tight

    
%% 3d Plot
%   view() tips:
%     AZ = -37.5, EL = 30 is the default 3-D view.
%     AZ = 0, EL = 90 is directly overhead and the default 2-D view.
%  
%     view(2) sets the default 2-D view, AZ = 0, EL = 90.
%     view(3) sets the default 3-D view, AZ = -37.5, EL = 30.
frames = 900;
% Creating view() vectors
% Azymuth Angle
lin_az = linspace(-37.5,0,frames);
% Elevation Angle
% lin_el = linspace(30,90,frames);
log_el = (logspace(0,2,frames) + 30) * (90/130);

% Gaussian curve function
gaussian = @(x,mean,dp) ((1./sqrt(2*pi*(dp.^2))) * exp(- ((x-mean).^2)./(2*dp.^2)));

% ------------------------------Id--------------------------------------- %
% Making the amplitude vector, in order to scan the time vs amp mesh
amp_vec = linspace(min(noisy_measure.id),max(noisy_measure.id), numel(time));

% Generate Density Probability function (z-axis) matrix
z_plt = zeros(numel(time));
for i = 1:numel(time)
    z_plt(i,:) = gaussian(amp_vec, kalman.id(i), dp.id(i));
end
% z_plt must be transposed!!!

figure()
id_plot = gcf;
% id_plot.Position = [640 0 1920 1080]; % Full-HD screen
id_plot.Position = [0 0 1360 768]; % HD screen
% Configuring Plot
    surf(time, amp_vec, z_plt', 'EdgeColor', 'none', 'FaceColor', 'interp',...
        'EdgeLighting','gouraud');
    daspect([1 33.3333 2433.3333])
    title('Kalman Estimative (Id)')
    xlabel('Time(s)')
    ylabel('Amplitude')
    zlabel('Density Probability Function')
    
% Create movie scene    
for i = 1:numel(lin_az)
    view(lin_az(i),log_el(i))
%     drawnow
    movievec(i) = getframe(id_plot);
end

mywriter = VideoWriter('Id_3d');
mywriter.FrameRate = 60;
mywriter.Quality = 100;
open(mywriter)
writeVideo(mywriter, movievec)
close(mywriter)
clear movievec

% ------------------------------Iq--------------------------------------- %
% Making the amplitude vector, in order to scan the time vs amp mesh
amp_vec = linspace(min(noisy_measure.iq),max(noisy_measure.iq), numel(time));

% Generate Density Probability function (z-axis) matrix
z_plt = zeros(numel(time));
for i = 1:numel(time)
    z_plt(i,:) = gaussian(amp_vec, kalman.iq(i), dp.iq(i));
end
% z_plt must be transposed!!!

figure()
iq_plot = gcf;
% iq_plot.Position = [640 0 1920 1080]; % Full-HD screen
iq_plot.Position = [0 0 1360 768]; % HD screen
% Configuring Plot
    surf(time, amp_vec, z_plt', 'EdgeColor', 'none', 'FaceColor', 'interp',...
        'EdgeLighting','gouraud');
    daspect([1 33.3333 2433.3333])
    title('Kalman Estimative (Iq)')
    xlabel('Time(s)')
    ylabel('Amplitude')
    zlabel('Density Probability Function')
    
% Create movie scene    
for i = 1:numel(lin_az)
    view(lin_az(i),log_el(i))
    drawnow
    movievec(i) = getframe(iq_plot);
end

mywriter = VideoWriter('Iq_3d');
mywriter.FrameRate = 60;
mywriter.Quality = 100;
open(mywriter)
writeVideo(mywriter, movievec)
close(mywriter)
clear movievec

% ------------------------------Omega----------------------------------- %
% Making the amplitude vector, in order to scan the time vs amp mesh
amp_vec = linspace(min(noisy_measure.omega ...
    ),max(noisy_measure.omega), numel(time));

% Generate Density Probability function (z-axis) matrix
z_plt = zeros(numel(time));
for i = 1:numel(time)
    z_plt(i,:) = gaussian(amp_vec, kalman.omega(i), dp.omega(i));
end
% z_plt must be transposed!!!

figure()
omega_plot = gcf;
% omega_plot.Position = [640 0 1920 1080]; % Full-HD screen
omega_plot.Position = [0 0 1360 768]; % HD screen
% Configuring Plot
    surf(time, amp_vec, z_plt', 'EdgeColor', 'none', 'FaceColor', 'interp',...
        'EdgeLighting','gouraud');
%     daspect([1 33.3333 2433.3333])
    title('Kalman Estimative (Omega)')
    xlabel('Time(s)')
    ylabel('Amplitude')
    zlabel('Density Probability Function')
    
% Create movie scene    
for i = 1:numel(lin_az)
    view(lin_az(i),log_el(i))
%     drawnow
    movievec(i) = getframe(omega_plot);
end

mywriter = VideoWriter('Omega_3d');
mywriter.FrameRate = 60;
mywriter.Quality = 100;
open(mywriter)
writeVideo(mywriter, movievec)
close(mywriter)
clear movievec

%     Generate GIF?
    % Seek help:
    % -view
    % -Representing Data as a Surface
    % -Low-Level Camera Properties