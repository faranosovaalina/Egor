%TEST_KIRCHHOFF_JET  Тесты kirchhoff_pse на PSE-решениях модельной струи.
%
%  1) St = 0.5, m = 0, M = 0.4, 0.9, 1.3, 1.5, 2: звуковое поле на полярной
%     дуге радиуса 200*R (R = 0.5 - радиус сопла) с центром на срезе сопла.
%  2) Ближнее поле: сравнение Кирхгофа с прямым PSE-давлением на линии
%     r = rLine недалеко от поверхности Кирхгофа (для больших M они должны
%     быть близки - PSE сам описывает излучение волн Маха).
%
% mat-файлы PSE создаются make_pse_mat.m, если их ещё нет (это самая долгая
% часть; в MATLAB - минуты).  Совместимо с MATLAB R2014a+ и Octave.

clear; close all;

Ms    = [0.4 0.9 1.3 1.5 2.0];
St    = 0.5;
m     = 0;
R     = 0.5;                       % радиус сопла, D = 1
gen   = [0.5 1.0; 14.0 4.0];       % образующая конуса (x0 r0; x1 r1)
theta = (15:5:165)';               % полярные углы, град от оси струи
Rarc  = 200*R;
obsFF = kirchhoff_arc(Rarc, theta, 0);
rLine = 5.0;                       % линия для сравнения в ближнем поле
% x >= 6: сюда волны Маха приходят с поверхности, а не из-за её начала x0
obsNF = [linspace(6, 13, 29)', rLine*ones(29, 1)];
kopts = struct('nSub', 4, 'taper', 0, 'verbose', false);

nM  = numel(Ms);
Dir = zeros(numel(theta), nM);     % |p|*Rarc на дуге
pK  = zeros(size(obsNF, 1), nM);   % Кирхгоф, ближнее поле
pP  = zeros(size(obsNF, 1), nM);   % PSE, ближнее поле
for i = 1:nM
    M = Ms(i);
    f = sprintf('pse_M%s_Tinf1p00_m%d_St_sweep.mat', strrep(sprintf('%.2f', M), '.', 'p'), m);
    if ~exist(f, 'file')
        % x0 = 0.5: у самого среза (x0 = 0.1) слой смешения модели тоньше
        % 0.015 D и мода К-Г не разрешается при Nr ~ 200
        make_pse_mat(M, St, m, struct('verbose', false, 'x0', 0.5, 'Nr', 160), f);
    end
    S  = load(f);
    iSt = find(abs(S.St_array - St) < 1e-9, 1);
    kopts.iSt = iSt;
    K  = kirchhoff_pse(S, gen, obsFF, kopts);
    Dir(:, i) = K.absp*Rarc;
    KN = kirchhoff_pse(S, gen, obsNF, kopts);
    pK(:, i) = KN.p;
    pP(:, i) = pse_pressure_at(S, iSt, obsNF);
    [~, ipk] = max(Dir(:, i));
    err = norm(pK(:, i) - pP(:, i))/norm(pP(:, i));
    fprintf(['M = %.1f: дальнее поле - максимум |p|*R при theta = %3d град; ' ...
        'ближнее поле (r = %.1f): отн. отличие Кирхгоф/PSE = %.3f\n'], ...
        M, theta(ipk), rLine, err);
end

%% ---- графики ------------------------------------------------------------
figure('Position', [100 100 1000 420]);
subplot(1, 2, 1);
semilogy(theta, Dir, 'LineWidth', 1.2); grid on;
xlabel('\theta, град'); ylabel('|p''| \cdot R_{arc}');
title(sprintf('Дуга R_{arc} = %g D, St = %.2f, m = %d', Rarc, St, m));
legend(arrayfun(@(v) sprintf('M = %.1f', v), Ms, 'UniformOutput', false), ...
    'Location', 'best');
subplot(1, 2, 2);
iM = nM;
plot(obsNF(:,1), real(pP(:, iM)), 'k-', obsNF(:,1), real(pK(:, iM)), 'r--', ...
     obsNF(:,1), abs(pP(:, iM)), 'k-', obsNF(:,1), abs(pK(:, iM)), 'r--', ...
     'LineWidth', 1.2);
grid on; xlabel('x/D');
title(sprintf('Ближнее поле r = %.1f, M = %.1f: PSE (черн.), Кирхгоф (красн.)', ...
    rLine, Ms(iM)));
