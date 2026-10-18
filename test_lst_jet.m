%TEST_LST_JET  Локальная задача устойчивости с деформацией контура (lst_jet) на
%модельной струе M = 0.9, St = 0.5, m = 0, и сравнение с alpha из PSE.
%
% Проверяется: (1) без деформации контура мода за нейтральным сечением
% "перескакивает" в непрерывный спектр (Im alpha -> 0, Re alpha -> const);
% (2) с деформацией она аналитически продолжается в Im alpha > 0;
% (3) результат за нейтральным сечением не зависит от ширины горба контура
% и сетки (пока мода не слилась с непрерывным спектром, dbranch > 0.03).
%
% Среднее течение - meanflow_model(..., smoothEps) с гладким переходом на
% границе ядра: у исходной модели профиль в этой точке не аналитичен, и
% продолжение в комплексную плоскость r невозможно (см. README).
% Нужен mat-файл PSE (его делает make_pse_mat.m).

clear; close all;
M = 0.9; St = 0.5; m = 0;
f = sprintf('pse_M%s_Tinf1p00_m%d_St_sweep.mat', strrep(sprintf('%.2f', M), '.', 'p'), m);
if ~exist(f, 'file')
    make_pse_mat(M, St, m, struct('verbose', false, 'x0', 0.5, 'Nr', 160), f);
end
S  = load(f);
mf = meanflow_model(M, 1, 1.4, 0.03);

base = struct('Nr', 140, 'rmax', 12, 'verbose', false, 'x', 0.5:0.25:6);
cfg   = {struct('saveFile', 'lst_M0p90_m0_St0p50.mat'), struct('widthScale', 1), ...
         struct('widthScale', 3), struct('Nr', 200), struct('deform', false)};
names = {'базовый', 'ширина x0.5', 'ширина x1.5', 'Nr = 200', 'БЕЗ деформации'};
A = cell(1, numel(cfg));
for c = 1:numel(cfg)
    o = base; fn = fieldnames(cfg{c});
    for j = 1:numel(fn), o.(fn{j}) = cfg{c}.(fn{j}); end
    out = lst_jet(mf, S, o);
    A{c} = out.alpha_all{1};
    if c == 1, out0 = out; end
    fprintf('%-16s x_neutral = %.3f, полюсов задето: %d\n', names{c}, out.x_neutral, ...
        sum(out.polesCrossed_all{1}));
end
x   = out0.x_all{1};
sel = x > out0.x_neutral + 0.3 & out0.dbranch_all{1} > 0.1;
fprintf('\nза нейтральным сечением, до слияния с непрерывным спектром (x = %.2f..%.2f):\n', ...
    min(x(sel)), max(x(sel)));
for c = 2:numel(cfg)
    fprintf('  max|alpha - alpha_базовый| (%-16s) = %.2e\n', names{c}, ...
        max(abs(A{c}(sel) - A{1}(sel))));
end

figure('Position', [100 100 900 380]);
subplot(1, 2, 1);
plot(x, real(A{1}), 'b-', x, real(A{end}), 'k:', x, real(out0.alphaPSE_all{1}), 'r--', 'LineWidth', 1.3);
grid on; xlabel('x/D'); ylabel('Re \alpha');
legend('LST, деформация', 'LST, без деформации', 'PSE', 'Location', 'best');
subplot(1, 2, 2);
plot(x, -imag(A{1}), 'b-', x, -imag(A{end}), 'k:', x, -imag(out0.alphaPSE_all{1}), 'r--', 'LineWidth', 1.3);
grid on; xlabel('x/D'); ylabel('-Im \alpha');
