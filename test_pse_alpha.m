%TEST_PSE_ALPHA  Варианты определения alpha в PSE после максимума пакета.
%
% Модельная струя (сглаженный излом, smoothEps = 0.03), St = 0.5, m = 0,
% M = 0.9 и 0.4. Варианты:
%   A - обычная нормировка по всем переменным (normVars = 1:5, как раньше);
%   B - нормировка по давлению (normVars = 'p');
%   C - alpha из LST (lst_jet) там, где оно действительно, далее - нормировка
%       по давлению.
% (Вариант rNorm = 2 тоже проверялся: он совпадает с A - причина второго
%  перехода Im alpha через ноль не во внешнем поле, а в конвективных
%  возмущениях ρ', T', u' внутри струи, см. README.)
% Проверяется: (1) сколько раз Im alpha меняет знак; (2) совпадает ли
% физическое поле p' около струи (r = 0.5, 1) во всех вариантах - выбор alpha
% меняет только разбиение p' = phat*exp(i*int alpha), поэтому различие
% показывает ошибку дискретизации; (3) "наблюдаемое" alpha_eff(x, r0).
% Совместимо с MATLAB R2014a+ и Octave. Время в Octave - ~10 мин на M.

clear; close all;
Ms = [0.9 0.4]; St = 0.5; m = 0;
rEff = [0.5 1.0];
base = struct('Nr', 140, 'rmax', 12, 'x0', 0.5, 'x1', 15, 'dx', 0.1, ...
              'rEff', rEff, 'verbose', false);
for iM = 1:numel(Ms)
    M  = Ms(iM);
    mf = meanflow_model(M, 1, 1.4, 0.03);
    o  = base; o.M = M;
    rA = pse_jet(mf, St, m, o);
    o.normVars = 'p';
    rB = pse_jet(mf, St, m, o);
    src = struct('x_all', {{rA.x}}, 'alpha_all', {{rA.alpha}}, 'St_array', St, ...
                 'm', m, 'M', M, 'Tinf', 1, 'gam', 1.4);
    L  = lst_jet(mf, src, struct('Nr', 140, 'rmax', 12, 'verbose', false));
    o.alphaPrescribed = L;
    rC = pse_jet(mf, St, m, o);

    R  = {rA, rB, rC}; names = {'A: норма 1:5', 'B: норма p', 'C: alpha LST'};
    x  = rA.x;
    fprintf('\nM = %.1f, St = %.2f:  LST: нейтр. x = %.2f, alpha_LST используется при x <= %.2f\n', ...
        M, St, L.x_neutral, max(x(rC.alphaSource == 1)));
    for v = 1:3
        a  = R{v}.alpha;
        sc = x(find(diff(sign(imag(a))) ~= 0) + 1);
        % поле около струи: сравнение с вариантом A
        pv = R{v}.p; pa = rA.p;
        jj = arrayfun(@(rr) find(abs(R{v}.r - rr) == min(abs(R{v}.r - rr)), 1), rEff);
        err = max(max(abs(pv(jj, :) - pa(jj, :))))/max(max(abs(pa(jj, :))));
        fprintf('  %-13s Im alpha меняет знак при x = %-18s  max|p''-p''_A|/max|p''_A| (r=0.5,1) = %.3f\n', ...
            names{v}, mat2str(round(100*sc)/100), err);
    end
    fprintf('     x  | Im a_A    Im a_B    Im a_C   | Im a_eff(r=1): A       B       C\n');
    for n = round(linspace(0.25*numel(x), numel(x), 8))
        fprintf('  %6.2f | %8.4f  %8.4f  %8.4f | %13.4f %7.4f %7.4f\n', x(n), ...
            imag(rA.alpha(n)), imag(rB.alpha(n)), imag(rC.alpha(n)), ...
            imag(rA.alphaEff(2, n)), imag(rB.alphaEff(2, n)), imag(rC.alphaEff(2, n)));
    end
    save(sprintf('test_pse_alpha_M%s.mat', strrep(sprintf('%.2f', M), '.', 'p')), ...
        'rA', 'rB', 'rC', 'L', '-v7');

    figure('Position', [100 100 1000 380]);
    subplot(1, 2, 1);
    plot(x, -imag(rA.alpha), 'k-', x, -imag(rB.alpha), 'b-', x, -imag(rC.alpha), 'r--', ...
         x, -imag(rB.alphaEff(2, :)), 'g:', 'LineWidth', 1.3);
    grid on; xlabel('x/D'); ylabel('-Im \alpha');
    legend('A', 'B', 'C', '\alpha_{eff}(r=1), B');
    title(sprintf('M = %.1f, St = %.2f', M, St));
    subplot(1, 2, 2);
    semilogy(x, abs(rA.p(jj(2), :)), 'k-', x, abs(rB.p(jj(2), :)), 'b-', ...
             x, abs(rC.p(jj(2), :)), 'r--', 'LineWidth', 1.3);
    grid on; xlabel('x/D'); ylabel('|p''(x, r = 1)|');
end
