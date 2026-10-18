function f = aaa_fit(Z, F, tol, mmax)
%AAA_FIT  Рациональная аппроксимация AAA (Nakatsukasa, Sete, Trefethen 2018)
%данных F(Z) на вещественных точках Z. Нужна для аналитического продолжения
%табличного профиля среднего течения в комплексную плоскость r.
%
%   f = aaa_fit(Z, F)            tol = 1e-10, mmax = 100
%   f = aaa_fit(Z, F, tol, mmax)
%   v = aaa_eval(f, z)           значения в комплексных точках z
%   p = aaa_poles(f)             полюса аппроксимации
%
% f: структура с узлами zj, значениями fj и весами wj барицентрической
% формы  R(z) = sum(wj*fj/(z-zj)) / sum(wj/(z-zj)),  f.err - ошибка на данных.

if nargin < 3 || isempty(tol),  tol  = 1e-10; end
if nargin < 4 || isempty(mmax), mmax = 100; end
Z = Z(:); F = F(:);
M = numel(Z);
J = (1:M).';
zj = zeros(0, 1); fj = zeros(0, 1);
C = zeros(M, 0);
R = mean(F)*ones(M, 1);
nF = max(abs(F));
if nF == 0, nF = 1; end
wj = 1; err = 0;
for m = 1:min(mmax, M-1)
    [~, k] = max(abs(F(J) - R(J)));
    jj = J(k);
    zj = [zj; Z(jj)]; fj = [fj; F(jj)];                    %#ok<AGROW>
    J(k) = [];
    C = [C, 1./(Z - Z(jj))];                               %#ok<AGROW>
    A = bsxfun(@times, bsxfun(@minus, F(J), fj.'), C(J, :));
    if size(A, 1) < m, break; end
    [~, ~, V] = svd(A, 0);
    wj = V(:, m);
    N = C*(wj.*fj); D = C*wj;
    R = F;
    R(J) = N(J)./D(J);
    err = max(abs(F - R));
    if err <= tol*nF, break; end
end
f.zj = zj; f.fj = fj; f.wj = wj; f.err = err/nF;
end
