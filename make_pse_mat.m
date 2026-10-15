function outFile = make_pse_mat(M, St, m, opts, outFile)
%MAKE_PSE_MAT  Расчёт PSE для модельной струи (та же модель среднего течения,
%что в example_new_mean_jet_to_mat.m) и сохранение mat-файла в том же
%формате (поля x_all, r_grid, St_array, P_all_St, alpha_all, M, Tinf, m, ...),
%который читает kirchhoff_pse.
%
%   outFile = make_pse_mat(M, St, m)
%   outFile = make_pse_mat(M, St, m, opts, outFile)
%
%   opts - параметры pse_jet (кроме M); дополнительно opts.Tinf (1),
%          opts.gamma (1.4). Значения по умолчанию - как в
%          example_new_mean_jet_to_mat.m.
%   Файл сохраняется в формате -v7 (читается и MATLAB, и Octave).

if nargin < 4 || isempty(opts), opts = struct(); end
def = struct('gamma', 1.4, 'Tinf', 1.0, 'Nr', 200, 'rmax', 12, 'rc', 0.5, ...
    'beta', 5, 'x0', 0.1, 'x1', 15, 'dx', 0.1, 'rOut', linspace(0, 6, 241)');
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i})), opts.(f{i}) = def.(f{i}); end
end
gam  = opts.gamma;
Tinf = opts.Tinf;
opts = rmfield(opts, 'Tinf');
opts.M = M;

mf  = meanflow_model(M, Tinf, gam);
res = pse_jet(mf, St, m, opts);

nSt = numel(res);
x_all = cell(1, nSt); P_all_St = cell(1, nSt); alpha_all = cell(1, nSt);
for k = 1:nSt
    x_all{k}     = res(k).x;
    P_all_St{k}  = res(k).p;
    alpha_all{k} = res(k).alpha;
end
r_grid   = res(1).r(:);
St_array = [res.St];
R  = 0.5;
Nr = numel(r_grid);

if nargin < 5 || isempty(outFile)
    tagM = strrep(sprintf('%.2f', M),    '.', 'p');
    tagT = strrep(sprintf('%.2f', Tinf), '.', 'p');
    outFile = sprintf('pse_M%s_Tinf%s_m%d_St_sweep.mat', tagM, tagT, m);
end
save(outFile, 'x_all', 'r_grid', 'St_array', 'P_all_St', 'alpha_all', ...
     'M', 'Tinf', 'm', 'gam', 'R', 'Nr', '-v7');
fprintf('MAT-файл сохранён: %s\n', outFile);
end

function mf = meanflow_model(M, Tinf, gam)
% модель среднего течения из example_new_mean_jet_to_mat.m (lambda = 0)
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
U = U/Uj0;
dUdx = zeros(size(U));
dUdx(2:end-1,:) = (U(3:end,:) - U(1:end-2,:))/(2*(x(2) - x(1)));
dUdx([1 end],:) = dUdx([2 end-1],:);
V = -cumtrapz(r, R.*dUdx, 2)./max(R, eps);  V(:,1) = 0;
T = Tinf + (1 - Tinf)*U + 0.5*(gam - 1)*M^2*U.*(1 - U);
mf = struct('x', x, 'r', r, 'U', U, 'V', V, 'T', T);
end
