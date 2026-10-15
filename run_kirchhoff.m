%RUN_KIRCHHOFF  Пример: звуковое поле по формуле Кирхгофа для всех St из
%mat-файла PSE (сохранённого example_new_mean_jet_to_mat.m / make_pse_mat.m).
% Совместимо с MATLAB R2014a+ и Octave.

clear; close all;

matFile = 'pse_M0p90_Tinf1p00_m0_St_sweep.mat';   % файл PSE
R       = 0.5;                                     % радиус сопла (D = 1)

% образующая конической поверхности Кирхгофа: (x0, r0) -> (x1, r1);
% должна лежать вне слоя смешения и внутри PSE-сетки (x в x_all, r <= max(r_grid))
gen = [0.5 1.0;
       14.0 4.0];

% точки наблюдения: полярная дуга радиуса 200R с центром на срезе сопла
theta = (10:2:170)';
obs   = kirchhoff_arc(200*R, theta, 0);
% (можно задать любые точки:  obs = [x(:), r(:)])

kopts = struct();
kopts.nSub  = 4;       % деление каждого шага PSE-сетки по x (1 - ровно PSE-сетка)
kopts.taper = 0;       % окно Тьюки у концов образующей, доля длины (напр. 0.05)
% kopts.Nphi = 256;    % число точек по азимуту (по умолчанию - авто)
kopts.saveFile = strrep(matFile, '.mat', '_kirchhoff.mat');

K = kirchhoff_pse(matFile, gen, obs, kopts);

% K.p [Nobs x nSt]: комплексное давление; K.absp, K.realp, K.imagp, K.phase
figure('Position', [100 100 1100 450]);
subplot(1, 2, 1);
semilogy(theta, K.absp, 'LineWidth', 1.2); grid on;
xlabel('\theta, град'); ylabel('|p''|');
title(sprintf('|p''| на дуге %g D, m = %d', 200*R, K.m));
legend(arrayfun(@(s) sprintf('St = %.2f', s), K.St, 'UniformOutput', false), ...
    'Location', 'best');
subplot(1, 2, 2);
plot(theta, unwrap(K.phase), 'LineWidth', 1.2); grid on;
xlabel('\theta, град'); ylabel('фаза p'', рад (unwrap)');
