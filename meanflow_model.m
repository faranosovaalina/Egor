function mf = meanflow_model(M, Tinf, gam, smoothEps)
%MEANFLOW_MODEL  Модель среднего течения круглой струи (та же, что в
%example_new_mean_jet_to_mat.m, lambda = 0): гауссов слой смешения вокруг
%потенциального ядра, Крокко-Буземан для температуры.
%   mf = meanflow_model(M, Tinf, gam)  ->  struct x, r, U, V, T
%   mf = meanflow_model(M, Tinf, gam, smoothEps)
%      smoothEps > 0 сглаживает излом профиля на границе ядра (r = rc) на
%      масштабе smoothEps, делая профиль аналитическим (нужно для проверки
%      деформации контура в lst_jet; по умолчанию 0 - исходная модель)
% модель среднего течения из example_new_mean_jet_to_mat.m (lambda = 0)
if nargin < 4 || isempty(smoothEps), smoothEps = 0; end
x = (0:0.05:20)';
r = [0:0.002:1.498, 1.5:0.01:3.99, 4:0.05:12]';
[R, X] = meshgrid(r, x);
lambda = 0; Uj0 = M;
Lc = 9.61*lambda + 0.775*(3*Uj0 + 5.13);
s1 = 0.65; s2 = 0.17; s3 = 0.41; L01 = 0.99; L02 = 1.25;
xi = X/Lc;
ds = zeros(size(X));
m1 = xi <= L01; m2 = xi > L01 & xi <= L02; m3 = xi > L02;
ds(m1) = s1*xi(m1);
ds(m2) = s1*L01 + s2*(xi(m2) - L01);
ds(m3) = s1*L01 + s2*(L02 - L01) + s3*(xi(m3) - L02);
ds = max(ds, 1e-6);
rc = max(0.5*(1 - X/Lc), 0);
Uj = Uj0*ones(size(X));
mw = X > 0.8*Lc;
Uj(mw) = Uj0*(1 - exp(-1./(0.7*(X(mw)/(0.8*Lc) - 1))));
U = zeros(size(R));
in = R <= rc; ou = ~in;
U(in) = Uj(in);
U(ou) = Uj(ou).*(lambda + (1 - lambda).*exp(-((R(ou) - rc(ou))./ds(ou)).^2));
if smoothEps > 0
    % гладкая версия: положительная часть (R - rc) через softplus
    z   = (R - rc)/smoothEps;
    pos = smoothEps*(max(z, 0) + log1p(exp(-abs(z))));
    U   = Uj.*(lambda + (1 - lambda).*exp(-(pos./ds).^2));
end
U = U/Uj0;
dUdx = zeros(size(U));
dUdx(2:end-1,:) = (U(3:end,:) - U(1:end-2,:))/(2*(x(2) - x(1)));
dUdx([1 end],:) = dUdx([2 end-1],:);
V = -cumtrapz(r, R.*dUdx, 2)./max(R, eps);  V(:,1) = 0;
T = Tinf + (1 - Tinf)*U + 0.5*(gam - 1)*M^2*U.*(1 - U);
mf = struct('x', x, 'r', r, 'U', U, 'V', V, 'T', T);
end
