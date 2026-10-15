%EXAMPLE_JET  Пример расчёта PSE для модельной круглой струи 
% "модельное среднее течение" на своё поле (x, r, U, V, rho/T),
% обезразмеренное по D_j, U_j, rho_j, T_j.

clear; close all;

%% ---- модельное среднее течение ----------------------------------------
M     = 0.4;          % акустическое число Маха струи U_j/c_j
gam   = 1.4;
Tinf  = 1.0;          % T_inf/T_j
x  = (0:0.05:20)';
r  = [0:0.002:1.498, 1.5:0.01:3.99, 4:0.05:12]';  % сгущение у слоя смешения
[R, X] = meshgrid(r, x);                       % массивы [Nx x Nr]
sp    = @(z) log(1 + exp(z));                  % гладкий max(z,0)
% xc    = 5;                                     % длина потенциального ядра
% theta = 0.02 + 0.03*X;                         % толщина слоя смешения
% rhalf = 0.5 + 0.08*sp(X - xc);                 % полуширина
% Uc    = 1./(1 + 0.12*sp(X - xc));              % скорость на оси
% F     = @(rr) 0.5*(1 + tanh((rhalf - rr)./(2*theta)));
% U     = Uc.*F(R)./F(0*R);
% % радиальная скорость из уравнения неразрывности (несжимаемое приближение)
% dUdx = zeros(size(U));
% dUdx(2:end-1,:) = (U(3:end,:) - U(1:end-2,:))/(2*(x(2)-x(1)));
% dUdx([1 end],:) = dUdx([2 end-1],:);
% V = -cumtrapz(r, R.*dUdx, 2)./max(R, eps);  V(:,1) = 0;
% % температура: соотношение Крокко-Буземана
% T = Tinf + (1 - Tinf)*U + 0.5*(gam - 1)*M^2*U.*(1 - U);
%% ---- модельное среднее течение (U-модель из find_modes_profiles_modelVc_Amp05.m) -------
% Uj0 = M (число Маха струи), в обезразмеривании PSE U_j = 1.
% lambda = Ucf/Uj0 = 0 (нет спутного потока).
lambda_mf = 0;
Uj0_mf    = M;

% Lc - длина начального участка (как в модели)
Lc_mf = 9.61*lambda_mf + 0.775*(3*Uj0_mf + 5.13);

% трёхкусочная модель толщины слоя смешения delta_s(x)
s1 = 0.65; s2 = 0.17; s3 = 0.41; L01 = 0.99; L02 = 1.25;
xi  = X / Lc_mf;
delta_s = zeros(size(X));
m1 = xi <= L01;
m2 = xi > L01 & xi <= L02;
m3 = xi > L02;
delta_s(m1) = s1*xi(m1);
delta_s(m2) = s1*L01 + s2*(xi(m2) - L01);
delta_s(m3) = s1*L01 + s2*(L02 - L01) + s3*(xi(m3) - L02);
delta_s = max(delta_s, 1e-6);

% полуширина ядра rc(x) = 0.5*(1 - x/Lc), обрезаем при rc < 0
rc = 0.5*(1 - X/Lc_mf);
rc = max(rc, 0);

% скорость на оси Uj(x) - модель Witze
Uj_x = Uj0_mf * ones(size(X));
mw = X > 0.8*Lc_mf;
Uj_x(mw) = Uj0_mf * (1 - exp(-1 ./ (0.7*(X(mw)/(0.8*Lc_mf) - 1))));

% продольная скорость U(x,r), гауссов профиль
inner = R <= rc;
outer = ~inner;
U = zeros(size(R));
U(inner) = Uj_x(inner);
U(outer) = Uj_x(outer) .* (lambda_mf + (1 - lambda_mf) .* ...
                           exp(-((R(outer) - rc(outer))./delta_s(outer)).^2));

% обезразмеривание на U_j (в PSE U_j = 1) - модели возвращают U/Uj0
U = U / Uj0_mf;

% --- радиальная скорость из уравнения неразрывности (как было) --------
dUdx = zeros(size(U));
dUdx(2:end-1,:) = (U(3:end,:) - U(1:end-2,:))/(2*(x(2)-x(1)));
dUdx([1 end],:) = dUdx([2 end-1],:);
V = -cumtrapz(r, R.*dUdx, 2)./max(R, eps);  V(:,1) = 0;

% --- температура: соотношение Крокко-Буземана (как было) --------------
T = Tinf + (1 - Tinf)*U + 0.5*(gam - 1)*M^2*U.*(1 - U);

mf = struct('x', x, 'r', r, 'U', U, 'V', V, 'T', T);

%mf = struct('x', x, 'r', r, 'U', U, 'V', V, 'T', T, 'rho', rho = 1./T);

%% ---- параметры PSE ------------------------------------------------------
St = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];      % числа Струхаля
m  = 0;                  % азимутальная мода
opts = struct();
opts.M     = M;
opts.gamma = gam;
opts.Nr    = 200;        % точек Чебышёва (для St > 0.5 у сопла нужно >= 140)
opts.rmax  = 12;         % внешняя граница
opts.rc    = 0.5;        % сгущение точек у слоя смешения
opts.beta  = 5;
opts.x0    = 0.1;        % начало маршевания (LST)
opts.x1    = 15;
opts.dx    = 0.1;
opts.rOut  = linspace(0, 6, 241)';   % сетка для вывода давления

res = pse_jet(mf, St, m, opts);

%% ---- сохранение в mat (ячейки на своих сетках) -------------------------
nSt = numel(res);

x_all     = cell(1, nSt);
P_all_St  = cell(1, nSt);
U_all     = cell(1, nSt);
Uj_1d     = cell(1, nSt);
A0_all    = cell(1, nSt);
Vc_all    = cell(1, nSt);
alpha_all = cell(1, nSt);

r_grid = res(1).r(:);          % общая (opts.rOut), проверяем ниже
Nr     = numel(r_grid);

for k = 1:nSt
    x_all{k}    = res(k).x;
    P_all_St{k} = res(k).p;

    % защита: r одинаков для всех St
    if numel(res(k).r) ~= Nr || max(abs(res(k).r(:) - r_grid)) > 1e-12
        warning('St = %.3g: r-сетка отличается от r_grid(1).', res(k).St);
    end

    % среднее U на (x_all{k}, r_grid)
    U_all{k} = interp_meanfield_U_local(mf, x_all{k}, r_grid);

    % Uj(x) из модели
    Uj_1d{k} = Uj_x_at(mf, x_all{k});

    % амплитуда на оси
    A0_all{k} = abs(res(k).p(1, :));

    % фазовая скорость и alpha
    Vc_all{k}    = res(k).omega ./ real(res(k).alpha);
    alpha_all{k} = res(k).alpha;
end

St_array = [res.St];

% --- общие скаляры ------------------------------------------------------
lambda = lambda_mf;
Uj0    = Uj0_mf;
Ucf    = 0;
Lc     = Lc_mf;
R      = 0.5;
R_max  = max(r_grid);
model   = struct();
vcmodel = struct();

% --- имя файла ----------------------------------------------------------
tagM = strrep(sprintf('%.2f', M),    '.', 'p');
tagT = strrep(sprintf('%.2f', Tinf), '.', 'p');
outFile = sprintf('pse_M%s_Tinf%s_m%d_St_sweep.mat', tagM, tagT, m);

save(outFile, ...
     'x_all', 'r_grid', 'St_array', ...
     'P_all_St', 'U_all', 'Uj_1d', 'A0_all', 'Vc_all', 'alpha_all', ...
     'lambda', 'Uj0', 'Ucf', 'Lc', 'R', 'R_max', ...
     'model', 'vcmodel', ...
     'M', 'Tinf', 'm', 'gam', 'Nr', ...
     '-v7.3');

fprintf('MAT-файл сохранён: %s\n', outFile);
fprintf('  nSt = %d, Nr = %d, Nx: ', nSt, Nr);
for k = 1:nSt, fprintf('%d ', numel(x_all{k})); end
fprintf('\n');

%% ---- графики ------------------------------------------------------------
% figure('Position', [100 100 1100 750]);
% for k = 1:numel(res)
%     subplot(numel(res), 2, 2*k-1);
%     pcolor(res(k).x, res(k).r, real(res(k).p)); shading interp;
%     axis equal tight; colorbar;
%     title(sprintf('Re p''(x,r), St = %.2f, m = %d', res(k).St, m));
%     xlabel('x/D'); ylabel('r/D');
%     subplot(numel(res), 2, 2*k);
%     yyaxis left;  semilogy(res(k).x, max(abs(res(k).p), [], 1)); ylabel('max_r |p''|');
%     yyaxis right; plot(res(k).x, -imag(res(k).alpha)); ylabel('-\alpha_i');
%     xlabel('x/D'); grid on;
% end
%% ---- графики ------------------------------------------------------------
figure('Position', [100 100 1100 750]);
for k = 1:numel(res)
    subplot(numel(res), 2, 2*k-1);
    pcolor(res(k).x, res(k).r, real(res(k).p)); shading interp;
    axis equal tight; colorbar;
    title(sprintf('Re p''(x,r), St = %.2f, m = %d', res(k).St, m));
    xlabel('x/D'); ylabel('r/D');
    
    subplot(numel(res), 2, 2*k);
    % --- две оси Y через plotyy (совместимо со старыми MATLAB) ---
    y1 = max(abs(res(k).p), [], 1);      % данные для левой оси (лог. шкала)
    y2 = -imag(res(k).alpha);            % данные для правой оси
    
    [hAx, hL1, hL2] = plotyy(res(k).x, y1, res(k).x, y2, 'semilogy', 'plot');
    ylabel(hAx(1), 'max_r |p''|');
    ylabel(hAx(2), '-\alpha_i');
    xlabel(hAx(1), 'x/D');
    grid(hAx(1), 'on'); grid(hAx(2), 'on');
    
    % Косметика: тонкие линии, чтобы совпадало с прежним видом
    set(hL1, 'LineWidth', 1.2);
    set(hL2, 'LineWidth', 1.2);
    set(hAx(1), 'YColor', 'k');   % можно поменять цвета осей
    set(hAx(2), 'YColor', 'k');
end
