function out = kirchhoff_pse(src, gen, obs, opts)
%KIRCHHOFF_PSE  Звуковое давление по формуле Кирхгофа (функция Грина
%свободного пространства) по полю давления PSE для одной азимутальной моды.
%
%   out = kirchhoff_pse(src, gen, obs)
%   out = kirchhoff_pse(src, gen, obs, opts)
%
% Совместимо с MATLAB R2014a+ и GNU Octave (>= 4.4).
%
% Поверхность Кирхгофа S - тело вращения вокруг оси струи x, заданное
% образующей в плоскости (x, r). Для конической поверхности образующая -
% отрезок  gen = [x0 r0; x1 r1].  Допускается и ломаная из нескольких
% вершин [x_1 r_1; x_2 r_2; ...] (например, коническая поверхность с
% торцами для замкнутой поверхности). Вершины перечисляются так, чтобы
% источник (струя) оставался "слева", т.е. внешняя нормаль n = (-t_r, t_x),
% где t - единичный касательный вектор образующей (для отрезка, идущего
% вниз по потоку, нормаль направлена наружу, в сторону больших r).
%
% Вне S среда неподвижна, звуковое поле удовлетворяет уравнению Гельмгольца
% с k = omega/c_inf, и (временной множитель exp(-i*omega*t))
%
%   p(x) = Int_S [ p(y) dG/dn_y - G dp/dn_y ] dS_y ,
%   G = exp(i*k*|x-y|) / (4*pi*|x-y|).
%
% На поверхности p(y) = p_m(x,r) exp(i*m*phi'); интеграл по phi' берётся
% численно (правило трапеций, спектрально точное для периодической функции),
% по образующей - правило трапеций. Результат дан в меридиональной плоскости
% phi = 0; для произвольного phi надо умножить на exp(i*m*phi).
%
% Разбиение образующей по x для каждого St берётся из исходной PSE-сетки
% x_all{k} (узлы PSE, попавшие в [x0, x1], плюс концы отрезка), каждый
% интервал дополнительно делится на opts.nSub частей. Значения p, dp/dx,
% dp/dr на образующей интерполируются с PSE-сетки. Если в файле есть
% alpha_all, интерполируется медленная огибающая p*exp(-i*int(alpha)dx) и
% отдельно фаза int(alpha)dx - это точнее прямой интерполяции осциллирующего p.
%
% ВХОД
%   src  - имя mat-файла (как сохраняет example_new_mean_jet_to_mat.m) или
%          структура с полями:
%            x_all     {1 x nSt}  x-сетки PSE (x/D)
%            r_grid    [Nr x 1]   r-сетка (r/D), общая для всех St
%            St_array  [1 x nSt]  числа Струхаля  St = f*D/U_j
%            P_all_St  {1 x nSt}  [Nr x Nx] комплексное давление PSE p(x,r)
%            alpha_all {1 x nSt}  alpha(x) (необязательно, но желательно)
%            M, Tinf, m           число Маха U_j/c_j, T_inf/T_j, азим. мода
%   gen  - образующая [x0 r0; x1 r1] (или ломаная, см. выше), в единицах D
%   obs  - [Nobs x 2] точки наблюдения (x, r) в плоскости phi = 0
%          (см. kirchhoff_arc - точки на полярной дуге)
%   opts - необязательные параметры:
%            .m       азимутальное число               (из файла, иначе 0)
%            .M       число Маха струи U_j/c_j          (из файла)
%            .Tinf    T_inf/T_j                         (из файла, иначе 1)
%            .cinf    скорость звука окр. среды в ед. U_j (sqrt(Tinf)/M)
%            .nSub    деление каждого интервала PSE-сетки по x   (4)
%            .Nphi    число точек по азимуту; [] - авто (по k и
%                     расстоянию до поверхности)               ([])
%            .taper   доля длины образующей у каждого открытого конца,
%                     на которой подынтегральное выражение плавно
%                     гасится окном Тьюки (0 - без окна)          (0)
%            .closed  поверхность замкнута (окно не применяется) (false)
%            .iSt     номера St из файла для расчёта      (все)
%            .verbose                                         (true)
%            .saveFile  имя mat-файла для сохранения результата ('')
%
% ВЫХОД  out:
%   .St    [1 x nSt],  .k [1 x nSt] = omega/c_inf,  .m, .cinf
%   .obs   [Nobs x 2]
%   .p     [Nobs x nSt] комплексное звуковое давление (ед. rho_j*U_j^2)
%   .absp, .realp, .imagp  - модуль, действительная, мнимая части
%   .phase [Nobs x nSt]    фаза angle(p), рад
%   .surf  {1 x nSt}       узлы поверхности: x, r, nx, nr, w, p, dpdn
%   .gen, .opts

if nargin < 4 || isempty(opts), opts = struct(); end
if ischar(src)
    S = load(src);
else
    S = src;
end

% --- данные PSE ---------------------------------------------------------
x_all = tocell(S.x_all);
P_all = tocell(S.P_all_St);
r     = S.r_grid(:);
St    = S.St_array(:).';
hasAlpha = isfield(S, 'alpha_all') && ~isempty(S.alpha_all);
if hasAlpha, alpha_all = tocell(S.alpha_all); end

m    = getopt(opts, 'm',    getfld(S, 'm', 0));
M    = getopt(opts, 'M',    getfld(S, 'M', NaN));
Tinf = getopt(opts, 'Tinf', getfld(S, 'Tinf', 1));
cinf = getopt(opts, 'cinf', sqrt(Tinf)/M);
nSub = getopt(opts, 'nSub', 4);
Nphi = getopt(opts, 'Nphi', []);
taper   = getopt(opts, 'taper', 0);
closed  = getopt(opts, 'closed', false);
iSt     = getopt(opts, 'iSt', 1:numel(St));
verbose = getopt(opts, 'verbose', true);
saveFile = getopt(opts, 'saveFile', '');
if ~isfinite(cinf) || cinf <= 0
    error('kirchhoff_pse:cinf', 'Не задана скорость звука: укажите opts.M или opts.cinf');
end
if size(gen, 2) ~= 2 || size(gen, 1) < 2
    error('kirchhoff_pse:gen', 'gen должен быть [x0 r0; x1 r1] (или ломаная N x 2)');
end
if size(obs, 2) ~= 2, error('kirchhoff_pse:obs', 'obs должен быть [Nobs x 2]: (x, r)'); end

nK   = numel(iSt);
Nobs = size(obs, 1);
out.St    = St(iSt);
out.k     = 2*pi*St(iSt)/cinf;
out.m     = m;
out.cinf  = cinf;
out.obs   = obs;
out.gen   = gen;
out.p     = zeros(Nobs, nK);
out.surf  = cell(1, nK);

for kk = 1:nK
    k  = iSt(kk);
    x  = x_all{k}(:).';
    P  = P_all{k};
    if size(P, 1) ~= numel(r) && size(P, 2) == numel(r), P = P.'; end
    omega = 2*pi*St(k);
    kw = omega/cinf;
    if hasAlpha
        alpha = alpha_all{k}(:).';
    else
        alpha = zeros(size(x));
    end
    sf = surface_field(x, r, P, alpha, gen, nSub);

    % плавное гашение у открытых концов образующей
    if taper > 0 && ~closed
        L  = sf.s(end);
        a  = taper*L;
        wt = ones(size(sf.s));
        i1 = sf.s < a;      wt(i1) = 0.5*(1 - cos(pi*sf.s(i1)/a));
        i2 = sf.s > L - a;  wt(i2) = 0.5*(1 - cos(pi*(L - sf.s(i2))/a));
        sf.w = sf.w.*wt;
    end

    if verbose
        fprintf('Кирхгоф: St = %.3g, m = %d, k = %.4g, узлов на образующей %d\n', ...
            St(k), m, kw, numel(sf.x));
    end
    out.p(:, kk) = kirchhoff_integral(sf, obs, kw, m, Nphi);
    out.surf{kk} = sf;
end

out.absp  = abs(out.p);
out.realp = real(out.p);
out.imagp = imag(out.p);
out.phase = angle(out.p);
out.opts  = struct('nSub', nSub, 'Nphi', Nphi, 'taper', taper, 'closed', closed);

if ~isempty(saveFile)
    St_K = out.St; k_K = out.k; obs_K = out.obs; p_K = out.p;
    absp = out.absp; realp = out.realp; imagp = out.imagp; phase = out.phase;
    gen_K = gen; m_K = m; cinf_K = cinf;
    save(saveFile, 'St_K', 'k_K', 'obs_K', 'p_K', 'absp', 'realp', 'imagp', ...
        'phase', 'gen_K', 'm_K', 'cinf_K', '-v7');
    if verbose, fprintf('Результат сохранён: %s\n', saveFile); end
end
end

% =========================================================================
%  Узлы на образующей и значения p, dp/dn в них
% =========================================================================
function sf = surface_field(x, r, P, alpha, gen, nSub)
% огибающая и фаза:  p = H*exp(i*Phi),  Phi = int alpha dx
Phi = cumtrapz(x, alpha);
H   = P.*repmat(exp(-1i*Phi), numel(r), 1);
Hx  = dmat(H, x, 2);
Hr  = dmat(H, r, 1);

xs = []; rs = []; nx = []; nr = []; w = []; s = []; s0 = 0;
for j = 1:size(gen, 1) - 1
    xa = gen(j, 1); ra = gen(j, 2); xb = gen(j+1, 1); rb = gen(j+1, 2);
    L  = hypot(xb - xa, rb - ra);
    if L == 0, continue; end
    % параметр t in [0,1]; базовые узлы - PSE-сетка по x (или по r для
    % вертикального отрезка), каждый интервал делится на nSub частей
    if abs(xb - xa) > 1e-12*L
        lo = min(xa, xb); hi = max(xa, xb);
        g  = x(x > lo + 1e-12 & x < hi - 1e-12);
        t  = (g - xa)/(xb - xa);
    else
        lo = min(ra, rb); hi = max(ra, rb);
        g  = r(r > lo + 1e-12 & r < hi - 1e-12).';
        t  = (g - ra)/(rb - ra);
    end
    t = unique([0, t(:).', 1]);
    t = refine(t, nSub);
    xt = xa + t*(xb - xa);
    rt = ra + t*(rb - ra);
    st = t*L;
    wt = zeros(size(t));                 % трапеции по длине дуги
    ds = diff(st);
    wt(1:end-1) = wt(1:end-1) + ds/2;
    wt(2:end)   = wt(2:end)   + ds/2;
    xs = [xs, xt]; rs = [rs, rt];                    %#ok<AGROW>
    nx = [nx, -(rb - ra)/L*ones(size(t))];           %#ok<AGROW>
    nr = [nr,  (xb - xa)/L*ones(size(t))];           %#ok<AGROW>
    w  = [w, wt]; s = [s, s0 + st];                  %#ok<AGROW>
    s0 = s0 + L;
end

if min(xs) < x(1) - 1e-9 || max(xs) > x(end) + 1e-9
    error('kirchhoff_pse:range', ['Образующая выходит за PSE-сетку по x: ' ...
        '[%.3g, %.3g] при x_PSE в [%.3g, %.3g]'], min(xs), max(xs), x(1), x(end));
end
if min(rs) < r(1) - 1e-9 || max(rs) > r(end) + 1e-9
    error('kirchhoff_pse:range', ['Образующая выходит за PSE-сетку по r: ' ...
        '[%.3g, %.3g] при r_PSE в [%.3g, %.3g]'], min(rs), max(rs), r(1), r(end));
end

Hs  = pse_interp_xr(x, r, H,  xs, rs);
Hxs = pse_interp_xr(x, r, Hx, xs, rs);
Hrs = pse_interp_xr(x, r, Hr, xs, rs);
Phs = cinterp1(x, Phi,   xs);
als = cinterp1(x, alpha, xs);
e   = exp(1i*Phs);

sf.x  = xs(:);  sf.r  = rs(:);
sf.nx = nx(:);  sf.nr = nr(:);
sf.w  = w(:);   sf.s  = s(:);
sf.p    = (Hs.*e).';
sf.p    = sf.p(:);
dpdx    = (Hxs + 1i*als.*Hs).*e;
dpdr    = Hrs.*e;
sf.dpdn = sf.nx.*dpdx(:) + sf.nr.*dpdr(:);
end

function t = refine(t, n)
if n <= 1, return; end
tt = zeros(1, (numel(t) - 1)*n + 1);
for i = 1:numel(t) - 1
    tt((i-1)*n + (1:n)) = t(i) + (0:n-1)/n*(t(i+1) - t(i));
end
tt(end) = t(end);
t = tt;
end

% =========================================================================
%  Интеграл Кирхгофа для набора точек наблюдения
% =========================================================================
function p = kirchhoff_integral(sf, obs, kw, m, NphiUser)
Nobs = size(obs, 1);
p = zeros(Nobs, 1);
xs = sf.x; rs = sf.r; nx = sf.nx; nr = sf.nr;
wr = sf.w.*rs;                                  % dS = r dphi ds
dsMax = max(diff(sf.s(sf.w > 0)));
if isempty(dsMax), dsMax = 0; end
warned = false;
for io = 1:Nobs
    X = obs(io, 1); rho = obs(io, 2);
    dmin = min(hypot(xs - X, rs - rho));
    if dmin < 1e-9
        error('kirchhoff_pse:obs', 'Точка наблюдения %d лежит на поверхности', io);
    end
    if ~warned && dmin < 2*dsMax
        warning('kirchhoff_pse:near', ['Точка наблюдения ближе к поверхности (%.3g), ' ...
            'чем 2 шага образующей (%.3g): увеличьте opts.nSub'], dmin, dsMax);
        warned = true;
    end
    if isempty(NphiUser)
        Nphi = max([64, ceil(4*(kw*max(rs) + abs(m)) + 32), ...
                    ceil(16*pi*max(rs)*(rho > 0)/dmin)]);
        Nphi = min(Nphi, 8192);
    else
        Nphi = NphiUser;
    end
    phi  = 2*pi*(0:Nphi-1)/Nphi;
    cph  = cos(phi);
    em   = exp(1i*m*phi)*(2*pi/Nphi);
    dx   = xs - X;
    A    = dx.^2 + rs.^2 + rho^2;
    D    = sqrt(max(bsxfun(@minus, A, (2*rho*rs)*cph), 0));
    G    = exp(1i*kw*D)./(4*pi*D);
    ydn  = bsxfun(@minus, dx.*nx + nr.*rs, (nr*rho)*cph);     % (y - x).n
    dGdn = (1i*kw - 1./D).*G.*ydn./D;
    I    = bsxfun(@times, dGdn, sf.p) - bsxfun(@times, G, sf.dpdn);
    p(io) = sum(wr.*(I*em.'));
end
end

% =========================================================================
%  Вспомогательные функции
% =========================================================================
function D = dmat(F, x, dim)
% производная 2-го порядка на неравномерной сетке вдоль измерения dim
if dim == 2, F = F.'; end
x = x(:); n = numel(x); D = zeros(size(F));
if n < 3
    D = repmat((F(end,:) - F(1,:))/(x(end) - x(1)), n, 1);
else
    h = diff(x);
    h1 = h(1:end-1); h2 = h(2:end);
    a = -h2./(h1.*(h1 + h2)); b = (h2 - h1)./(h1.*h2); c = h1./(h2.*(h1 + h2));
    D(2:n-1,:) = bsxfun(@times, a, F(1:n-2,:)) + bsxfun(@times, b, F(2:n-1,:)) + ...
                 bsxfun(@times, c, F(3:n,:));
    h1 = h(1); h2 = h(2);
    D(1,:) = -(2*h1+h2)/(h1*(h1+h2))*F(1,:) + (h1+h2)/(h1*h2)*F(2,:) - h1/(h2*(h1+h2))*F(3,:);
    h1 = h(n-1); h2 = h(n-2);
    D(n,:) = (2*h1+h2)/(h1*(h1+h2))*F(n,:) - (h1+h2)/(h1*h2)*F(n-1,:) + h1/(h2*(h1+h2))*F(n-2,:);
end
if dim == 2, D = D.'; end
end

function v = cinterp1(x, f, xq)
v = interp1(x, real(f), xq, 'spline') + 1i*interp1(x, imag(f), xq, 'spline');
end

function c = tocell(a)
if iscell(a), c = a; else c = {a}; end
end

function v = getopt(s, name, def)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else v = def; end
end

function v = getfld(s, name, def)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else v = def; end
end
