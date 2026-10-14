%EXAMPLE_JET  Пример расчёта PSE для модельной круглой струи M = 0.9.
% Замените блок "модельное среднее течение" на своё поле (x, r, U, V, rho/T),
% обезразмеренное по D_j, U_j, rho_j, T_j.

clear; close all;

%% ---- модельное среднее течение ----------------------------------------
M     = 0.9;          % акустическое число Маха струи U_j/c_j
gam   = 1.4;
Tinf  = 1.0;          % T_inf/T_j
x  = (0:0.05:20)';
r  = [0:0.002:1.498, 1.5:0.01:3.99, 4:0.05:12]';  % сгущение у слоя смешения
[R, X] = meshgrid(r, x);                       % массивы [Nx x Nr]
sp    = @(z) log(1 + exp(z));                  % гладкий max(z,0)
xc    = 6;                                     % длина потенциального ядра
theta = 0.02 + 0.03*X;                         % толщина слоя смешения
rhalf = 0.5 + 0.08*sp(X - xc);                 % полуширина
Uc    = 1./(1 + 0.12*sp(X - xc));              % скорость на оси
F     = @(rr) 0.5*(1 + tanh((rhalf - rr)./(2*theta)));
U     = Uc.*F(R)./F(0*R);
% радиальная скорость из уравнения неразрывности (несжимаемое приближение)
dUdx = zeros(size(U));
dUdx(2:end-1,:) = (U(3:end,:) - U(1:end-2,:))/(2*(x(2)-x(1)));
dUdx([1 end],:) = dUdx([2 end-1],:);
V = -cumtrapz(r, R.*dUdx, 2)./max(R, eps);  V(:,1) = 0;
% температура: соотношение Крокко-Буземана
T = Tinf + (1 - Tinf)*U + 0.5*(gam - 1)*M^2*U.*(1 - U);

mf = struct('x', x, 'r', r, 'U', U, 'V', V, 'T', T);   % rho = 1./T

%% ---- параметры PSE ------------------------------------------------------
St = [0.3 0.5 0.7];      % числа Струхаля
m  = 0;                  % азимутальная мода
opts = struct();
opts.M     = M;
opts.gamma = gam;
opts.Nr    = 140;        % точек Чебышёва (для St > 0.5 у сопла нужно >= 140)
opts.rmax  = 12;         % внешняя граница
opts.rc    = 0.5;        % сгущение точек у слоя смешения
opts.beta  = 5;
opts.x0    = 0.2;        % начало маршевания (LST)
opts.x1    = 15;
opts.dx    = 0.1;
opts.rOut  = linspace(0, 6, 241)';   % сетка для вывода давления

res = pse_jet(mf, St, m, opts);

%% ---- графики ------------------------------------------------------------
figure('Position', [100 100 1100 750]);
for k = 1:numel(res)
    subplot(numel(res), 2, 2*k-1);
    pcolor(res(k).x, res(k).r, real(res(k).p)); shading interp;
    axis equal tight; colorbar;
    title(sprintf('Re p''(x,r), St = %.2f, m = %d', res(k).St, m));
    xlabel('x/D'); ylabel('r/D');
    subplot(numel(res), 2, 2*k);
    yyaxis left;  semilogy(res(k).x, max(abs(res(k).p), [], 1)); ylabel('max_r |p''|');
    yyaxis right; plot(res(k).x, -imag(res(k).alpha)); ylabel('-\alpha_i');
    xlabel('x/D'); grid on;
end
