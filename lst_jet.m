function out = lst_jet(mf, src, opts)
%LST_JET  Локальная (плоскопараллельная) задача устойчивости круглой струи в
%заданных сечениях x с отслеживанием моды Кельвина-Гельмгольца и деформацией
%радиального контура интегрирования вокруг критического слоя.
%
%   out = lst_jet(mf, src)
%   out = lst_jet(mf, src, opts)
%
% Совместимо с MATLAB R2014a+ и GNU Octave.
%
% В каждом сечении x среднее течение "замораживается": V = 0, d/dx = 0,
% профили U(r), T(r), rho = 1/T. Решается та же система линеаризованных
% уравнений Эйлера, что в pse_jet (переменные [rho,u_x,u_r,u_theta,T]),
% с возмущениями exp(i(alpha*x + m*theta - omega*t)): обобщённая задача на
% собственные значения по alpha
%       (-i*omega*A + i*m*D + E + B d/dr) q = alpha * (-i*C) q .
% Размерности и нормировки - как в pse_jet (D_j, U_j, rho_j, T_j, p = rho*T/(gamma*M^2)).
%
% КРИТИЧЕСКИЙ СЛОЙ. Уравнения сингулярны там, где omega - alpha*U(r) = 0, т.е.
% при U(r_c) = omega/alpha. Для неустойчивой моды (Im alpha < 0) Im(omega/alpha) > 0,
% и, т.к. в слое смешения U убывает с r, точка r_c лежит в нижней полуплоскости
% комплексного r (контур - вещественная ось - проходит выше неё). Когда
% Im alpha -> 0 (нейтральное сечение), r_c подходит к вещественной оси и
% переходит в верхнюю полуплоскость; если контур оставить вещественным, то он
% "перескочит" через особенность, и найденное решение перестанет быть
% аналитическим продолжением решения из неустойчивой области (мода
% превратилась бы в моду непрерывного спектра). Чтобы этого не случилось,
% контур деформируется:
%       r(s) = s + i*h(s),   s in [0, rmax] (сетка Чебышёва по s),
%       h(s) = amp * B(s),
% (B - гладкая "шапка" с плато вокруг Re r_c и компактным носителем) так,
% что в точке s = Re r_c контур остаётся по ту же сторону от r_c, что и
% вначале, на расстоянии не менее d:
%       sigma*(h(Re r_c) - Im r_c) >= d,   sigma = +1 (контур выше r_c).
% Положение r_c(x) находится из U(r_c) = omega/alpha по текущему alpha;
% амплитуда h гладко (softplus) включается по мере приближения r_c к
% контуру, так что вдали от нейтрального сечения контур вещественный.
%
% АНАЛИТИЧЕСКОЕ ПРОДОЛЖЕНИЕ ПРОФИЛЯ. Среднее течение задано таблицей на
% вещественных r, поэтому для комплексных r профили U(r) и T(r) в каждом
% сечении продолжаются рациональной аппроксимацией AAA (aaa_fit.m). Если
% профиль имеет излом/разрыв (как гауссов профиль с ядром в модели), AAA
% ставит полюса рядом с вещественной осью; их положение контролируется:
% out.polesCrossed показывает, сколько полюсов оказалось между вещественной
% осью и деформированным контуром (должно быть 0, иначе продолжение
% недостоверно). Можно задать собственные аналитические профили через
% opts.profileFun (см. ниже).
%
% СОПРОВОЖДЕНИЕ МОДЫ. Первое сечение: собственное значение, ближайшее к alpha
% PSE (из src) в этом сечении (или к opts.alphaGuess). Далее - метод Ньютона
% для пары (alpha, q) с начальным приближением из предыдущего сечения; если
% Ньютон не сошёлся или alpha скакнул более чем на opts.maxJump, то
% берётся ближайшее к ожидаемому значение из полного спектра.
%
% ВХОД
%   mf   - среднее течение (как в pse_jet): x [Nx], r [Nr], U [Nx x Nr],
%          T и/или rho [Nx x Nr]  (V не используется)
%   src  - mat-файл (или структура) расчёта PSE: x_all, alpha_all, St_array,
%          m, M, Tinf, gam. Определяет сечения (x_all каждого St), начальные
%          приближения alpha, M, m. Можно не давать (src = []), тогда
%          нужны opts.St, opts.x, opts.M, opts.m, opts.alphaGuess.
%   opts - необязательные параметры:
%            .iSt        номера St из src                        (все)
%            .St, .x     явные St и сечения x (вместо src)
%            .M, .m, .Tinf, .gamma  (из src)
%            .Nr         число точек Чебышёва по r                 (140)
%            .rmax       внешняя граница                           (12)
%            .mapping    'sinh' | 'malik'                          ('sinh')
%            .rc, .beta  сгущение 'sinh' (центр, интенсивность)    (0.5, 5)
%            .rmap       параметр 'malik'                          (1.0)
%            .deform     деформировать контур                      (true)
%            .margin     d / (толщина слоя смешения delta)         (0.3)
%            .widthScale ширина w / delta                          (2)
%            .side       sigma = +1/-1 (по умолч. определяется в 1-м сечении)
%            .contourIter макс. число пересчётов контура в сечении (4)
%            .aaaTol, .aaaPts, .aaaMax  параметры AAA            (1e-10, 500, 120)
%            .profileFun @(x, r) -> [U, T]  аналитические профили при комплексном r
%                        (если задан, AAA не используется)
%            .maxJump    допустимый скачок alpha между сечениями   (1.0)
%            .alphaGuess начальное приближение alpha в 1-м сечении
%            .newtonTol  точность Ньютона                          (1e-10)
%            .storeQ     сохранять полные собственные векторы      (false)
%            .saveFile   имя mat-файла для сохранения              ('')
%            .verbose                                              (true)
%
% ВЫХОД out (элементы 1 x nSt, ячейки по St):
%   .St, .m, .M, .Tinf, .gamma
%   .x_all{k}      [1 x Nx]  сечения
%   .alpha_all{k}  [1 x Nx]  alpha(x) локальной задачи (аналитически продолженная)
%   .alphaPSE_all{k}         alpha PSE в тех же сечениях (если есть src)
%   .P_all_St{k}   [Nr x Nx] собственная функция давления на контуре
%                            (нормировка p = 1 в точке максимума |p|)
%   .r_contour_all{k} [Nr x Nx] комплексный контур r(s)
%   .r_grid        [Nr x 1]  вещественный параметр s контура
%   .rcrit_all{k}  [1 x Nx]  критический слой r_c (комплексный)
%   .amp_all{k}    [1 x Nx]  амплитуда деформации контура h
%   .delta_all{k}  [1 x Nx]  оценка толщины слоя смешения
%   .dbranch_all{k} [1 x Nx] |alpha - omega/U_max|/(omega/U_max): близость к точке ветвления
%                            непрерывного спектра (при < 0.03 мода сливается с ним)
%   .continuum_all{k} [1 x Nx] true - Im alpha = 0 с машинной точностью: найдено
%                            значение непрерывного спектра, а не мода (недостоверно)
%   .width_all{k}  [1 x Nx]  ширина w горба контура (сужается, чтобы не задевать полюса профиля)
%   .polesCrossed_all{k} [1 x Nx], .newton_all{k} [1 x Nx] (1 - сошёлся Ньютон)
%   .x_neutral(k)  x, где Im alpha меняет знак (NaN, если не меняет)
%   .sigma         сторона обхода критического слоя

if nargin < 3 || isempty(opts), opts = struct(); end
if nargin < 2, src = []; end
if ischar(src), src = load(src); end
hasSrc = isstruct(src) && isfield(src, 'x_all');

m    = getopt(opts, 'm',    getfld(src, 'm', 0));
M    = getopt(opts, 'M',    getfld(src, 'M', NaN));
Tinf = getopt(opts, 'Tinf', getfld(src, 'Tinf', 1));
gam  = getopt(opts, 'gamma', getfld(src, 'gam', 1.4));
if ~isfinite(M), error('lst_jet:M', 'Не задано число Маха: opts.M или src.M'); end
kap = 1/(gam*M^2);

g.Nr = getopt(opts, 'Nr', 140);
g.rmax = getopt(opts, 'rmax', 12);
g.mapping = getopt(opts, 'mapping', 'sinh');
g.rc = getopt(opts, 'rc', 0.5);
g.beta = getopt(opts, 'beta', 5);
g.rmap = getopt(opts, 'rmap', 1.0);
deform  = getopt(opts, 'deform', true);
marginR = getopt(opts, 'margin', 0.3);
widthS  = getopt(opts, 'widthScale', 2);
side    = getopt(opts, 'side', 0);
cIter   = getopt(opts, 'contourIter', 4);
aaaTol  = getopt(opts, 'aaaTol', 1e-10);
aaaPts  = getopt(opts, 'aaaPts', 500);
aaaMax  = getopt(opts, 'aaaMax', 120);
profFun = getopt(opts, 'profileFun', []);
maxJump = getopt(opts, 'maxJump', 1.0);
ntol    = getopt(opts, 'newtonTol', 1e-10);
storeQ  = getopt(opts, 'storeQ', false);
saveFile = getopt(opts, 'saveFile', '');
verbose = getopt(opts, 'verbose', true);

% --- среднее течение -------------------------------------------------
mf.x = mf.x(:); mf.r = mf.r(:);
hasT = isfield(mf, 'T') && ~isempty(mf.T);
if ~hasT, mf.T = 1./mf.rho; end
if g.rmax > mf.r(end)
    warning('lst_jet:rmax', 'rmax уменьшен до mf.r(end) = %.3g', mf.r(end));
    g.rmax = mf.r(end);
end
idr = find(mf.r <= g.rmax + 1e-12);
step = max(1, ceil(numel(idr)/aaaPts));
idr = idr(1:step:end);
if idr(end) ~= find(mf.r <= g.rmax + 1e-12, 1, 'last')
    idr(end+1) = find(mf.r <= g.rmax + 1e-12, 1, 'last');
end
rd = mf.r(idr);

% --- St и сечения ------------------------------------------------------
if hasSrc
    Sall = src.St_array(:).';
    iSt  = getopt(opts, 'iSt', 1:numel(Sall));
    Stv  = Sall(iSt);
else
    Stv = opts.St(:).'; iSt = 1:numel(Stv);
end
nSt = numel(Stv);

% --- сетка по r --------------------------------------------------------
[sr, De, drde] = radial_grid(g.Nr, g);
N = numel(sr);
[bcIdx, bcRows] = boundary_rows(N, m, De);   % строки ГУ (De(1,:) не масштабируем - см. ниже)

out.St = Stv; out.m = m; out.M = M; out.Tinf = Tinf; out.gamma = gam;
out.r_grid = sr;
out.x_all = cell(1, nSt); out.alpha_all = cell(1, nSt); out.alphaPSE_all = cell(1, nSt);
out.P_all_St = cell(1, nSt); out.r_contour_all = cell(1, nSt);
out.rcrit_all = cell(1, nSt); out.amp_all = cell(1, nSt); out.delta_all = cell(1, nSt);
out.polesCrossed_all = cell(1, nSt); out.width_all = cell(1, nSt); out.dbranch_all = cell(1, nSt); out.continuum_all = cell(1, nSt); out.newton_all = cell(1, nSt);
out.Q_all = cell(1, nSt);
out.x_neutral = nan(1, nSt); out.sigma = zeros(1, nSt);

for kk = 1:nSt
    St = Stv(kk); omega = 2*pi*St;
    if hasSrc
        k = iSt(kk);
        xs = getopt(opts, 'x', src.x_all{k}(:).');
        xp = src.x_all{k}(:).';
        if isfield(src, 'alpha_all') && ~isempty(src.alpha_all)
            ap = src.alpha_all{k}(:).';
            aPSE = interp1(xp, real(ap), xs, 'linear', 'extrap') + ...
                   1i*interp1(xp, imag(ap), xs, 'linear', 'extrap');
        else
            aPSE = nan(size(xs));
        end
    else
        xs = opts.x(:).'; aPSE = nan(size(xs));
    end
    nx = numel(xs);
    if verbose, fprintf('LST: St = %.4g, m = %d, %d сечений\n', St, m, nx); end

    alpha = zeros(1, nx); P = zeros(N, nx); RC = zeros(N, nx);
    rcr = nan(1, nx) + 1i*nan(1, nx); amp = zeros(1, nx); dlt = zeros(1, nx);
    npole = zeros(1, nx); nwt = zeros(1, nx); wused = zeros(1, nx); dbr = zeros(1, nx);
    if storeQ, Qs = zeros(5*N, nx); end
    sigma = side; qprev = []; rcprev = NaN;

    for n = 1:nx
        xq = min(max(xs(n), mf.x(1)), mf.x(end));
        Ud = interp1(mf.x, mf.U, xq); Td = interp1(mf.x, mf.T, xq);
        Ud = Ud(idr).'; Td = Td(idr).';
        dU = gradient_nu(Ud, rd);
        delta = (max(Ud) - min(Ud))/max(abs(dU));
        if isempty(profFun)
            fU = aaa_fit(rd, Ud, aaaTol, aaaMax);
            fT = aaa_fit(rd, Td, aaaTol, aaaMax);
            pol = [aaa_poles(fU); aaa_poles(fT)];
        else
            fU = []; fT = []; pol = [];
        end
        if n == 1
            if isfield(opts, 'alphaGuess') && ~isempty(opts.alphaGuess)
                a0 = opts.alphaGuess(min(kk, numel(opts.alphaGuess)));
            elseif isfinite(aPSE(1))
                a0 = aPSE(1);
            else
                error('lst_jet:guess', 'Нужен opts.alphaGuess или alpha_all в src');
            end
        else
            if n >= 3
                a0 = alpha(n-1) + (alpha(n-1) - alpha(n-2))*(xs(n) - xs(n-1))/(xs(n-1) - xs(n-2));
                if abs(a0 - alpha(n-1)) > 0.5*maxJump, a0 = alpha(n-1); end
            else
                a0 = alpha(n-1);
            end
        end

        a = a0; ampPrev = inf; used = 0;
        for it = 1:cIter
            % --- контур по текущему alpha ------------------------------
            rcrit = NaN;
            if deform
                [rcrit, sigma] = crit_layer(fU, profFun, xq, rd, Ud, omega/a, rcprev, sigma);
            end
            d = marginR*delta; w = widthS*delta;
            if isfinite(rcrit)
                sc = min(max(real(rcrit), 0), g.rmax);
                [~, isc] = min(abs(sr - sc));
                wMin = 4*max(abs(diff(sr(max(isc-1,1):min(isc+1,N)))));
                % ширина сужается, пока контур не пересекает полюса профиля
                for sh = 1:25
                    [ampn, sc] = contour_amp(rcrit, sigma, d, w, g.rmax);
                    [r, drdeta, Dr, hh] = make_contour(sr, De, drde, ampn, sc, w, g.rmax);
                    if isempty(pol) || count_crossed(pol, sr, hh) == 0 || w <= wMin, break; end
                    w = max(0.8*w, wMin);
                end
            else
                ampn = 0; sc = 0.5*g.rmax;
                [r, drdeta, Dr, hh] = make_contour(sr, De, drde, ampn, sc, w, g.rmax);
            end
            [L, R, rho, T] = build_system(r, Dr, fU, fT, profFun, xq, omega, m, kap, gam);
            L(bcIdx, :) = bcRows(Dr);
            R(bcIdx, :) = 0;
            % --- собственное значение -----------------------------------
            [a, q, used] = solve_eig(L, R, a, qprev, ntol, maxJump);
            if abs(ampn - ampPrev) < 1e-3*d, break; end
            ampPrev = ampn;
        end
        if isfinite(rcrit), rcprev = rcrit; end
        % --- нормировка по давлению ---------------------------------
        p = kap*(rho.*q(4*N+1:5*N) + T.*q(1:N));
        [~, imx] = max(abs(p));
        q = q/p(imx); p = p/p(imx);
        qprev = q;
        alpha(n) = a; P(:, n) = p; RC(:, n) = r;
        rcr(n) = rcrit; amp(n) = ampn; dlt(n) = delta; nwt(n) = used;
        if ~isempty(pol), npole(n) = count_crossed(pol, sr, hh); end
        wused(n) = w;
        dbr(n) = abs(a - omega/max(Ud))/(omega/max(Ud));
        if storeQ, Qs(:, n) = q; end
        if verbose && (n == 1 || mod(n, 10) == 0 || n == nx)
            fprintf('  x = %6.3f: alpha = %8.4f %+8.4fi  (PSE %8.4f %+8.4fi)  r_c = %6.3f %+7.3fi  h = %6.3f\n', ...
                xs(n), real(a), imag(a), real(aPSE(n)), imag(aPSE(n)), real(rcrit), imag(rcrit), ampn);
        end
    end

    ineg = find(diff(sign(imag(alpha))) ~= 0, 1);
    if ~isempty(ineg)
        ai = imag(alpha(ineg:ineg+1));
        out.x_neutral(kk) = xs(ineg) - ai(1)*(xs(ineg+1) - xs(ineg))/(ai(2) - ai(1));
    end
    if verbose
        fprintf('  нейтральное сечение (Im alpha = 0): x = %.3f\n', out.x_neutral(kk));
        imrg = find(dbr < 0.03 & (1:nx) > 1, 1);
        if ~isempty(imrg)
            fprintf(['  x >= %.3f: alpha в пределах 3%% от omega/U_max (точка ветвления непрерывного\n' ...
                     '  спектра ядра) - мода сливается с непрерывным спектром, далее результат\n' ...
                     '  зависит от контура и физического смысла не имеет\n'], xs(imrg));
        end
        icn = find(abs(imag(alpha)) < 1e-6*abs(alpha));
        if ~isempty(icn)
            warning('lst_jet:continuum', ['в %d сечениях (x = %s) Im alpha = 0 с машинной ' ...
                'точностью: решение попало на непрерывный спектр, а не на моду К-Г ' ...
                '(continuum_all = true)'], numel(icn), mat2str(xs(icn), 3));
        end
        if max(npole) > 0
            warning('lst_jet:poles', ['в %d сечениях контур пересекает полюса аналитического продолжения ' ...
                'профиля (излом/разрыв профиля?) - результат там недостоверен'], sum(npole > 0));
        end
    end
    out.x_all{kk} = xs; out.alpha_all{kk} = alpha; out.alphaPSE_all{kk} = aPSE;
    out.P_all_St{kk} = P; out.r_contour_all{kk} = RC; out.rcrit_all{kk} = rcr;
    out.amp_all{kk} = amp; out.delta_all{kk} = dlt;
    out.continuum_all{kk} = abs(imag(alpha)) < 1e-6*abs(alpha); out.dbranch_all{kk} = dbr; out.width_all{kk} = wused; out.polesCrossed_all{kk} = npole; out.newton_all{kk} = nwt; out.sigma(kk) = sigma;
    if storeQ, out.Q_all{kk} = Qs; end
end

if ~isempty(saveFile)
    x_all = out.x_all; alpha_all = out.alpha_all; alphaPSE_all = out.alphaPSE_all; %#ok<NASGU>
    P_all_St = out.P_all_St; r_contour_all = out.r_contour_all; r_grid = out.r_grid; %#ok<NASGU>
    rcrit_all = out.rcrit_all; amp_all = out.amp_all; delta_all = out.delta_all; %#ok<NASGU>
    polesCrossed_all = out.polesCrossed_all; width_all = out.width_all; dbranch_all = out.dbranch_all; continuum_all = out.continuum_all; x_neutral = out.x_neutral; %#ok<NASGU>
    sigma = out.sigma; St_array = out.St; Nr = N; gam = out.gamma; %#ok<NASGU>
    save(saveFile, 'x_all', 'alpha_all', 'alphaPSE_all', 'P_all_St', ...
        'r_contour_all', 'r_grid', 'rcrit_all', 'amp_all', 'delta_all', ...
        'polesCrossed_all', 'width_all', 'dbranch_all', 'continuum_all', 'x_neutral', 'sigma', 'St_array', 'm', 'M', ...
        'Tinf', 'gam', 'Nr', '-v7');
    if verbose, fprintf('LST: результат сохранён в %s\n', saveFile); end
end
end

% =========================================================================
%  Критический слой и контур
% =========================================================================
function [rc, sigma] = crit_layer(fU, profFun, xq, rd, Ud, c, rguess, sigma)
% решает U(r) = c в комплексной плоскости r (Ньютон)
cr = real(c);
if cr >= max(Ud) || cr <= min(Ud)
    rc = NaN; return
end
ufun = @(z) prof_U(fU, profFun, xq, z);
r0 = rguess;
if ~isfinite(r0)
    k = find(Ud >= cr, 1, 'last');
    if isempty(k) || k >= numel(rd), rc = NaN; return; end
    r0 = rd(k) + (Ud(k) - cr)/(Ud(k) - Ud(k+1))*(rd(k+1) - rd(k));
    r0 = r0 + 0i;
end
hh = 1e-6;
r = r0;
for it = 1:60
    f  = ufun(r) - c;
    fp = (ufun(r + hh) - ufun(r - hh))/(2*hh);
    dr = -f/fp;
    if abs(dr) > 0.2, dr = 0.2*dr/abs(dr); end
    r = r + dr;
    if abs(dr) < 1e-12, break; end
end
rc = r;
if sigma == 0
    if abs(imag(rc)) > 1e-8, sigma = sign(-imag(rc)); else sigma = 1; end
end
end

function nc = count_crossed(pol, sr, hh)
% число полюсов, лежащих между вещественной осью и контуром r = s + i*h(s)
hint = interp1(sr, hh, real(pol), 'linear', NaN);
nc = sum(isfinite(hint) & (imag(pol).*hint > 0) & (abs(imag(pol)) < abs(hint)));
end

function u = prof_U(fU, profFun, xq, z)
if isempty(profFun), u = aaa_eval(fU, z);
else
    [u, ~] = profFun(xq, z);
end
end

function [amp, sc] = contour_amp(rcrit, sigma, d, ~, ~)
% амплитуда деформации в точке sc = Re r_c: контур должен остаться на расстоянии
% не менее d от r_c по ту же сторону (sigma = +1: выше), включается гладко (softplus)
sc = max(real(rcrit), 0);
t = sigma*imag(rcrit) + d;
e = d/2;
if t/e > 30, sp = t; else sp = e*log(1 + exp(t/e)); end
amp = sigma*sp;
end

function [r, drdeta, Dr, hh] = make_contour(sr, De, drde, amp, sc, w, rmax)
% h(s) = amp * B(s), B - гладкая (C-бесконечная) "шапка" с плато: B = 1 на
% [sc - wp, sc + wp], плавно спадает до 0 на длине w (но не шире, чем
% позволяют ось и внешняя граница), вне носителя h = 0 точно.
sc = min(sc, rmax);
wp = min(0.5*w, 0.25*max(min(sc, rmax - sc), 1e-3));
lL = min(w, 0.5*sc);            % длина спада слева
lR = min(w, 0.5*(rmax - sc));   % длина спада справа
B  = ones(size(sr));
if lL > 0, B = B.*smooth_step((sr - (sc - wp - lL))/lL); end
if lR > 0, B = B.*(1 - smooth_step((sr - (sc + wp))/lR)); end
hh = amp*B;
hp = (De*hh)./drde;
r  = sr + 1i*hh;
drdeta = drde.*(1 + 1i*hp);
Dr = bsxfun(@rdivide, De, drdeta);
end

function y = smooth_step(t)
% C-бесконечная ступенька: 0 при t <= 0, 1 при t >= 1
y = zeros(size(t));
y(t >= 1) = 1;
i = t > 0 & t < 1;
ti = t(i);
y(i) = 1./(1 + exp(1./ti - 1./(1 - ti)));
end

% =========================================================================
%  Система линеаризованных уравнений для параллельного течения
% =========================================================================
function [L, R, rho, T] = build_system(r, Dr, fU, fT, profFun, xq, omega, m, kap, gam)
N = numel(r);
if isempty(profFun)
    U = aaa_eval(fU, r); T = aaa_eval(fT, r);
else
    [U, T] = profFun(xq, r);
end
U = U(:); T = T(:);
rho = 1./T;
Ur = Dr*U; Tr = Dr*T; rhor = Dr*rho;
ri = zeros(N, 1); ri(2:N) = 1./r(2:N);       % строка оси заменяется ГУ
g1 = gam - 1;

c = cell(5); c{1,1} = 1; c{2,2} = rho; c{3,3} = rho; c{4,4} = rho; c{5,5} = 1;
A = assemble(c, N);
c = cell(5);
c{1,1} = U;        c{1,2} = rho;
c{2,1} = kap*T;    c{2,2} = rho.*U;  c{2,5} = kap*rho;
c{3,3} = rho.*U;   c{4,4} = rho.*U;
c{5,2} = g1*T;     c{5,5} = U;
C = assemble(c, N);
c = cell(5);
c{1,3} = rho;
c{3,1} = kap*T;    c{3,5} = kap*rho;
c{5,3} = g1*T;
B = assemble(c, N);
c = cell(5);
c{1,4} = rho.*ri;
c{4,1} = kap*T.*ri; c{4,5} = kap*rho.*ri;
c{5,4} = g1*T.*ri;
D = assemble(c, N);
c = cell(5);
c{1,3} = rhor + rho.*ri;
c{2,3} = rho.*Ur;
c{3,1} = kap*Tr;   c{3,5} = kap*rhor;
c{5,3} = Tr + g1*T.*ri;
E = assemble(c, N);

L = -1i*omega*A + 1i*m*D + E + B*kron(eye(5), Dr);
R = -1i*C;
end

function Mx = assemble(c, N)
Mx = zeros(5*N);
for i = 1:5
    for j = 1:5
        if ~isempty(c{i,j})
            d = c{i,j};
            if isscalar(d), d = d*ones(N, 1); end
            Mx((i-1)*N + (1:N), (j-1)*N + (1:N)) = diag(d(:));
        end
    end
end
end

function [idx, rowsFun] = boundary_rows(N, m, ~)
% индексы заменяемых строк и функция, строящая их по Dr (нужна Dr(1,:))
id = @(v, j) (v-1)*N + j;
idxF = id((1:5)', N);
switch abs(m)
    case 0
        idxA = [id(1,1); id(2,1); id(3,1); id(4,1); id(5,1)];
    case 1
        idxA = [id(1,1); id(2,1); id(3,1); id(4,1); id(5,1)];
    otherwise
        idxA = id((1:5)', 1);
end
idx = [idxA; idxF];
rowsFun = @(Dr) bc_rows(N, m, Dr, idxA, idxF);
end

function rows = bc_rows(N, m, Dr, idxA, idxF)
I  = eye(5*N);
id = @(v, j) (v-1)*N + j;
blk = @(v, row) [zeros(1,(v-1)*N), row, zeros(1,(5-v)*N)];
switch abs(m)
    case 0
        ra = [blk(1, Dr(1,:)); blk(2, Dr(1,:)); I(id(3,1),:); I(id(4,1),:); blk(5, Dr(1,:))];
    case 1
        ra = [I(id(1,1),:); I(id(2,1),:); blk(3, Dr(1,:)); ...
              I(id(3,1),:) + 1i*m*I(id(4,1),:); I(id(5,1),:)];
    otherwise
        ra = I(idxA, :);
end
rows = [ra; I(idxF, :)];
end

% =========================================================================
%  Собственное значение: Ньютон с продолжением, запасной вариант - eig
% =========================================================================
function [a, q, used] = solve_eig(L, R, a0, qprev, ntol, maxJump)
used = 0; a = a0; q = [];
n = size(L, 1);
if ~isempty(qprev)
    e = qprev; qq = qprev/(e'*qprev); aa = a0; ok = false;
    for it = 1:12
        Rq = R*qq;
        res = [(L - aa*R)*qq; e'*qq - 1];
        J = [L - aa*R, -Rq; e', 0];
        dd = J\(-res);
        qq = qq + dd(1:n); aa = aa + dd(end);
        if abs(dd(end)) < ntol*max(1, abs(aa)) && norm(dd(1:n)) < 1e-8*norm(qq)
            ok = true; break
        end
        if ~isfinite(aa), break; end
    end
    if ok && abs(aa - a0) < maxJump
        a = aa; q = qq; used = 1; return
    end
end
[V, Dg] = eig(L, R);
ev = diag(Dg);
good = isfinite(ev) & abs(ev) < 1e3;
cand = find(good);
[~, j] = min(abs(ev(cand) - a0));
j = cand(j);
a = ev(j); q = V(:, j);
end

% =========================================================================
%  Сетка по r (как в pse_jet) и вспомогательные функции
% =========================================================================
function [r, De, drde] = radial_grid(N, g)
[D0, eta] = cheb(N-1);
eta = flipud(eta); De = rot90(D0, 2);
switch lower(g.mapping)
    case 'malik'
        a = g.rmap*g.rmax/(g.rmax - 2*g.rmap);
        b = 1 + 2*a/g.rmax;
        r = a*(1 + eta)./(b - eta);
        drde = a*(1 + b)./(b - eta).^2;
    case 'sinh'
        rc = g.rc; be = g.beta; q = rc/g.rmax;
        A  = log((1 + (exp(be) - 1)*q)/(1 + (exp(-be) - 1)*q))/(2*be);
        xi = (1 + eta)/2;
        r  = rc*(1 + sinh(be*(xi - A))/sinh(be*A));
        drde = 0.5*rc*be*cosh(be*(xi - A))/sinh(be*A);
    otherwise
        error('lst_jet:grid', 'Неизвестное отображение %s', g.mapping);
end
r(1) = 0; r(end) = g.rmax;
end

function [D, x] = cheb(N)
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x, 1, N+1);
dX = X - X';
D = (c*(1./c)')./(dX + eye(N+1));
D = D - diag(sum(D, 2));
end

function d = gradient_nu(f, x)
f = f(:); x = x(:);
d = zeros(size(f));
d(2:end-1) = (f(3:end) - f(1:end-2))./(x(3:end) - x(1:end-2));
d(1) = (f(2) - f(1))/(x(2) - x(1));
d(end) = (f(end) - f(end-1))/(x(end) - x(end-1));
end

function v = getopt(s, name, def)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else v = def; end
end

function v = getfld(s, name, def)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else v = def; end
end
