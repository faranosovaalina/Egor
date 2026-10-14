function res = pse_jet(mf, St, m, opts)
%PSE_JET  Линейные параболизованные уравнения устойчивости (PSE) для круглой
%струи, по постановке Gudmundsson (PhD thesis, Caltech, 2009;
%Gudmundsson & Colonius, JFM 2011).
%
%   res = pse_jet(mf, St, m, opts)
%
% Решаются линеаризованные уравнения Эйлера для сжимаемого газа в
% цилиндрических координатах (x, r, theta) относительно возмущений
%   q' = [rho', u_x', u_r', u_theta', T']
% с анзацем PSE
%   q'(x,r,theta,t) = qhat(x,r) * exp( i*( int_x0^x alpha(xi) dxi + m*theta - omega*t ) ),
% где qhat - медленно меняющаяся по x функция формы, alpha(x) - комплексное
% волновое число. Неоднозначность разбиения "амплитуда/фаза" снимается
% условием нормировки
%   int_0^rmax qhat^H * d(qhat)/dx * r dr = 0 .
% Начальное условие - собственная функция Кельвина-Гельмгольца локальной
% (параллельной) задачи устойчивости в сечении x0. Маршевание по x -
% неявная схема Эйлера; по r - коллокация Чебышёва с отображением.
%
% Обезразмеривание (как у Gudmundsson): длины - диаметром сопла D_j,
% скорости - U_j, плотность - rho_j, температура - T_j, давление - rho_j*U_j^2.
% Уравнение состояния:  p = rho*T/(gamma*M^2),  M = U_j/c_j  (акуст. число Маха струи).
% Число Струхаля St = f*D_j/U_j,  omega = 2*pi*St.
%
% ВХОД
%   mf   - структура со средним течением на прямоугольной сетке:
%            mf.x   [Nx x 1]   осевые координаты (возрастают)
%            mf.r   [Nr x 1]   радиальные координаты (возрастают, желательно с r=0)
%            mf.U   [Nx x Nr]  осевая скорость
%            mf.V   [Nx x Nr]  радиальная скорость (необязательно, по умолч. 0)
%            mf.rho [Nx x Nr]  плотность    } достаточно задать одно из двух:
%            mf.T   [Nx x Nr]  температура  } при p = const  rho*T = 1
%   St   - вектор чисел Струхаля
%   m    - азимутальное волновое число (0, 1, 2, ...)
%   opts - структура параметров (все поля необязательны):
%            .M         число Маха струи U_j/c_j                    (0.9)
%            .gamma     показатель адиабаты                         (1.4)
%            .Nr        число точек Чебышёва по r                   (120)
%            .rmax      внешняя граница расчётной области          (15)
%            .mapping   отображение по r: 'sinh' или 'malik'          ('sinh')
%            .rc, .beta  'sinh': центр сгущения и его интенсивность    (0.5, 5)
%            .rmap      'malik': половина точек в [0, rmap]          (1.0)
%            .x0, .x1   начало и конец маршевания            (mf.x(1), mf.x(end))
%            .dx        шаг маршевания                              (0.1)
%            .dxAuto    увеличивать dx до dxSafety/Re(alpha0)       (true)
%            .dxSafety  коэффициент в ограничении dx > 1/Re(alpha)  (1.1)
%            .dropDpdx  отбросить d(phat)/dx (Li & Malik), снимает
%                       ограничение на шаг                          (false)
%            .tol       точность итераций по alpha                  (1e-8)
%            .maxIter   макс. число итераций по alpha на шаге       (30)
%            .normVars  переменные в условии нормировки, подмнож. 1:5 (1:5)
%            .alphaGuess  вектор начальных приближений alpha(x0) для
%                       каждого St (NaN -> автоматический выбор моды K-H) ([])
%            .cphRange  диапазон фазовой скорости моды K-H, доли max(U) ([0.2 2])
%            .minGrowth мин. -Im(alpha)/Re(alpha) моды К-Г при авто-выборе (0.02)
%            .lstTol    допуск совпадения собств. значений на двух сетках
%                       (отсев ложных мод Чебышёва), отн.             (1e-2)
%            .amp       амплитуда max_r|phat| в сечении x0           (1)
%            .interp    метод интерполяции среднего течения         ('pchip')
%            .storeQ    сохранять все компоненты qhat               (false)
%            .rOut      радиальная сетка для вывода давления        ([] - сетка Чебышёва)
%            .verbose                                               (true)
%
% ВЫХОД  res(k), k = 1..numel(St):
%   .St, .omega, .m
%   .x      [1 x Nx]   сечения маршевания
%   .r      [Nr x 1]   радиальная сетка (или opts.rOut)
%   .alpha  [1 x Nx]   комплексное волновое число alpha(x)
%   .phase  [1 x Nx]   int_x0^x alpha dx
%   .phat   [Nr x Nx]  функция формы давления
%   .p      [Nr x Nx]  комплексная амплитуда давления p'(x,r) =
%                      phat*exp(i*int alpha dx); физическое поле
%                      Re{ p * exp(i*m*theta - i*omega*t) }
%   .q      [5*Nr x Nx] (если storeQ) все функции формы
%   .converged [1 x Nx] признак сходимости итераций по alpha
%   .alphaLST           собственное значение в x0

if nargin < 4, opts = struct(); end
opts = set_defaults(opts, mf, numel(St));
mf   = prepare_meanflow(mf);

gam = opts.gamma;
kap = 1/(gam*opts.M^2);

% --- радиальная сетка и операторы ---------------------------------------
[r, Dr, wq] = radial_grid(opts.Nr, opts);
N  = numel(r);
DR = kron(eye(5), Dr);
% веса квадратуры для  int f r dr
W  = repmat(wq(:).*r(:), 5, 1);
wmask = zeros(5*N,1);
for v = opts.normVars(:)', wmask((v-1)*N+(1:N)) = 1; end
W = W.*wmask;
ip = @(a,b) sum(W.*conj(a).*b);
[bcIdx, bcRows] = boundary_rows(N, m, Dr);

for k = 1:numel(St)
    omega = 2*pi*St(k);
    if opts.verbose
        fprintf('PSE: St = %.4g, m = %d\n', St(k), m);
    end

    % --- 1. локальная задача устойчивости в x0 --------------------------
    s0 = meanflow_at(mf, opts.x0, r, opts.interp);
    s0p = make_parallel(s0);
    [alpha0, q0] = lst_kh(s0p, omega, m, opts, kap, gam, opts.alphaGuess(k));
    if opts.verbose
        fprintf('  LST:  alpha(x0) = %.5f %+.5fi  (c_ph = %.3f)\n', ...
            real(alpha0), imag(alpha0), omega/real(alpha0));
    end

    % --- 2. шаг маршевания ----------------------------------------------
    dx = opts.dx;
    if ~opts.dropDpdx
        dxmin = opts.dxSafety/abs(real(alpha0));
        if dx < dxmin
            if opts.dxAuto
                if opts.verbose
                    fprintf('  dx увеличен с %.4g до %.4g (условие dx > 1/Re(alpha))\n', dx, dxmin);
                end
                dx = dxmin;
            else
                warning('pse_jet:dx', ['dx = %.3g < 1/Re(alpha) = %.3g: возможна ' ...
                    'численная неустойчивость PSE.'], dx, 1/abs(real(alpha0)));
            end
        end
    end
    nst = max(2, ceil((opts.x1 - opts.x0)/dx) + 1);
    x   = linspace(opts.x0, opts.x1, nst);
    dx  = x(2) - x(1);

    % --- 3. нормировка начального условия по давлению ------------------
    ph0 = kap*(s0.rho.*q0(4*N+1:5*N) + s0.T.*q0(1:N));
    [~, imx] = max(abs(ph0));
    q0 = q0*(opts.amp/ph0(imx));

    % --- 4. маршевание --------------------------------------------------
    Q     = zeros(5*N, nst);  Q(:,1) = q0;
    alpha = zeros(1, nst);    alpha(1) = alpha0;
    conv  = true(1, nst);
    PH    = zeros(N, nst);    PH(:,1) = kap*(s0.rho.*q0(4*N+1:5*N) + s0.T.*q0(1:N));

    for n = 2:nst
        s = meanflow_at(mf, x(n), r, opts.interp);
        [A, B, C, D, E] = lee_matrices(s, r, kap, gam);
        L0 = -1i*omega*A + 1i*m*D + E + B*DR;
        Cx = C;                         % коэффициент при d(qhat)/dx
        if opts.dropDpdx
            Cx(N+1:2*N, 1:N)       = 0; % d(phat)/dx в осевом импульсе
            Cx(N+1:2*N, 4*N+1:5*N) = 0;
        end
        qp  = Q(:, n-1);
        rhs = Cx*qp/dx;
        rhs(bcIdx) = 0;
        % невязка условия нормировки  F(a) = <q, dq/dx>/<q, q> = 0:
        % первый шаг - стандартная поправка a - i*F, далее метод секущих
        % (простая итерация сходится медленно/расходится при малом dx)
        a = alpha(n-1);
        ok = false; aOld = []; FOld = [];
        for it = 1:opts.maxIter
            Mat = L0 + 1i*a*C + Cx/dx;
            Mat(bcIdx, :) = bcRows;
            qn = Mat\rhs;
            F  = ip(qn, (qn - qp)/dx)/ip(qn, qn);
            if isempty(aOld) || F == FOld
                da = -1i*F;
            else
                da = -F*(a - aOld)/(F - FOld);
            end
            aOld = a; FOld = F;
            a  = a + da;
            if abs(da) < opts.tol*max(1, abs(a)), ok = true; break; end
        end
        if ~ok && opts.verbose
            fprintf('  x = %.3f: итерации по alpha не сошлись (|dalpha| = %.2e)\n', x(n), abs(da));
        end
        conv(n)  = ok;
        alpha(n) = aOld;                % alpha, отвечающее найденному qn
        Q(:, n)  = qn;
        PH(:, n) = kap*(s.rho.*qn(4*N+1:5*N) + s.T.*qn(1:N));
    end

    phase = cumtrapz(x, alpha);
    P = PH .* repmat(exp(1i*phase), N, 1);

    rr = r;
    if ~isempty(opts.rOut)
        rr = opts.rOut(:);
        PH = interp_r(r, PH, rr);
        P  = interp_r(r, P,  rr);
    end

    res(k).St    = St(k);                 %#ok<AGROW>
    res(k).omega = omega;
    res(k).m     = m;
    res(k).x     = x;
    res(k).r     = rr;
    res(k).alpha = alpha;
    res(k).alphaLST = alpha0;
    res(k).phase = phase;
    res(k).phat  = PH;
    res(k).p     = P;
    res(k).converged = conv;
    if opts.storeQ
        res(k).q = Q; res(k).rq = r;
    end
    if opts.verbose
        [~, ipk] = max(max(abs(P), [], 1));
        fprintf('  макс. |p| при x = %.2f; не сошлось шагов: %d из %d\n', ...
            x(ipk), sum(~conv), nst-1);
    end
end
end

% =========================================================================
%  Линеаризованные уравнения Эйлера:
%    -i*omega*A q + i*alpha*C q + i*m*D q + E q + B dq/dr + C dq/dx = 0
%  (D уже содержит множитель 1/r).  q = [rho; u; v; w; T], блоки по N.
% =========================================================================
function [A, B, C, D, E] = lee_matrices(s, r, kap, gam)
N = numel(r);
ri = zeros(N,1); ri(r > 0) = 1./r(r > 0);   % на оси строки заменяются ГУ
g1 = gam - 1;
U = s.U; V = s.V; rho = s.rho; T = s.T;
divU = s.Ux + s.Vr + V.*ri;

c = cell(5); c{1,1} = 1; c{2,2} = rho; c{3,3} = rho; c{4,4} = rho; c{5,5} = 1;
A = assemble(c, N);

c = cell(5);                       % d/dx
c{1,1} = U;        c{1,2} = rho;
c{2,1} = kap*T;    c{2,2} = rho.*U;  c{2,5} = kap*rho;
c{3,3} = rho.*U;
c{4,4} = rho.*U;
c{5,2} = g1*T;     c{5,5} = U;
C = assemble(c, N);

c = cell(5);                       % d/dr
c{1,1} = V;        c{1,3} = rho;
c{2,2} = rho.*V;
c{3,1} = kap*T;    c{3,3} = rho.*V;  c{3,5} = kap*rho;
c{4,4} = rho.*V;
c{5,3} = g1*T;     c{5,5} = V;
B = assemble(c, N);

c = cell(5);                       % (1/r) d/dtheta
c{1,4} = rho.*ri;
c{4,1} = kap*T.*ri; c{4,5} = kap*rho.*ri;
c{5,4} = g1*T.*ri;
D = assemble(c, N);

c = cell(5);                       % члены без производных от q
c{1,1} = divU;   c{1,2} = s.rx;   c{1,3} = s.rr + rho.*ri;
c{2,1} = U.*s.Ux + V.*s.Ur + kap*s.Tx;
c{2,2} = rho.*s.Ux;  c{2,3} = rho.*s.Ur;  c{2,5} = kap*s.rx;
c{3,1} = U.*s.Vx + V.*s.Vr + kap*s.Tr;
c{3,2} = rho.*s.Vx;  c{3,3} = rho.*s.Vr;  c{3,5} = kap*s.rr;
c{4,4} = rho.*V.*ri;
c{5,2} = s.Tx;   c{5,3} = s.Tr + g1*T.*ri;  c{5,5} = g1*divU;
E = assemble(c, N);
end

function M = assemble(c, N)
M = zeros(5*N);
for i = 1:5
    for j = 1:5
        if ~isempty(c{i,j})
            d = c{i,j};
            if isscalar(d), d = d*ones(N,1); end
            M((i-1)*N+(1:N), (j-1)*N+(1:N)) = diag(d(:));
        end
    end
end
end

% =========================================================================
%  Граничные условия: ось (r = 0, индекс 1) и внешняя граница (индекс N)
% =========================================================================
function [idx, rows] = boundary_rows(N, m, Dr)
I = eye(5*N);
id = @(v, j) (v-1)*N + j;
blkrow = @(v, row) [zeros(1,(v-1)*N), row, zeros(1,(5-v)*N)];
idx = []; rows = [];
% ось
switch abs(m)
    case 0      % v = w = 0, d(rho,u,T)/dr = 0
        idx  = [id(1,1); id(2,1); id(3,1); id(4,1); id(5,1)];
        rows = [blkrow(1, Dr(1,:)); blkrow(2, Dr(1,:)); I(id(3,1),:); ...
                I(id(4,1),:); blkrow(5, Dr(1,:))];
    case 1      % rho = u = T = 0, dv/dr = 0, v + i*m*w = 0
        idx  = [id(1,1); id(2,1); id(3,1); id(4,1); id(5,1)];
        rows = [I(id(1,1),:); I(id(2,1),:); blkrow(3, Dr(1,:)); ...
                I(id(3,1),:) + 1i*m*I(id(4,1),:); I(id(5,1),:)];
    otherwise   % все компоненты равны нулю
        idx  = id((1:5)', 1);
        rows = I(idx, :);
end
% внешняя граница: однородные условия Дирихле
idxF = id((1:5)', N);
idx  = [idx; idxF];
rows = [rows; I(idxF, :)];
end

% =========================================================================
%  Локальная (параллельная) задача устойчивости: обобщённая задача на
%  собственные значения  L0 q = alpha*(-i C) q.  Выбирается мода
%  Кельвина-Гельмгольца; ложные моды отсеиваются сравнением с решением
%  на сетке другого размера.
% =========================================================================
function [alpha, q] = lst_kh(s, omega, m, opts, kap, gam, guess)
[ev1, q1] = lst_eig(s.src, opts.Nr, omega, m, opts, kap, gam, true);
ev2       = lst_eig(s.src, round(0.8*opts.Nr), omega, m, opts, kap, gam, false);

% физичные кандидаты: конечные, с умеренным |Im alpha|/Re alpha; при
% автоматическом выборе - заметно растущие (отсев почти нейтральных мод
% непрерывного спектра alpha ~ omega/U)
ok = isfinite(ev1) & abs(ev1) < 1e3 & abs(imag(ev1)) < 3*abs(real(ev1));
if isnan(guess)
    ok = ok & -imag(ev1) > opts.minGrowth*abs(real(ev1));
end
match = false(size(ev1));
for j = find(ok)'
    match(j) = min(abs(ev2 - ev1(j))) < opts.lstTol*abs(ev1(j));
end

if ~isnan(guess)
    cand = find(ok & match);
    if isempty(cand), cand = find(ok); end
    [~, j] = min(abs(ev1(cand) - guess));
    j = cand(j);
else
    Umax = max(s.U);
    cph  = omega./real(ev1);
    band = ok & real(ev1) > 0 & cph > opts.cphRange(1)*Umax & ...
           cph < opts.cphRange(2)*Umax & imag(ev1) < 0;
    cand = find(band & match);
    if isempty(cand)
        cand = find(band);
        if isempty(cand)
            error('pse_jet:lst', ['Не найдена неустойчивая мода К-Г при St = %.3g. ' ...
                'Задайте opts.alphaGuess или измените cphRange/x0.'], omega/(2*pi));
        end
        warning('pse_jet:lst', 'Мода К-Г не подтверждена на второй сетке; проверьте Nr.');
    end
    [~, j] = min(imag(ev1(cand)));
    j = cand(j);
end
alpha = ev1(j);
q = q1(:, j);
end

function [ev, V] = lst_eig(src, Nr, omega, m, opts, kap, gam, wantVec)
[r, Dr] = radial_grid(Nr, opts);
s = make_parallel(meanflow_at(src.mf, src.x, r, opts.interp));
[A, B, C, D, E] = lee_matrices(s, r, kap, gam);
L = -1i*omega*A + 1i*m*D + E + B*kron(eye(5), Dr);
R = -1i*C;
[idx, rows] = boundary_rows(Nr, m, Dr);
L(idx, :) = rows;
R(idx, :) = 0;
if wantVec
    [V, Dg] = eig(L, R);
    ev = diag(Dg);
else
    ev = eig(L, R);
    V = [];
end
end

function s = make_parallel(s)
s.V(:) = 0; s.Ux(:) = 0; s.Vx(:) = 0; s.Vr(:) = 0; s.rx(:) = 0; s.Tx(:) = 0;
end

% =========================================================================
%  Радиальная сетка: точки Гаусса-Лобатто-Чебышёва eta in [-1,1], r(-1)=0,
%  r(1)=rmax.  Два варианта отображения:
%   'sinh'  - сгущение точек около слоя смешения r = rc (интенсивность beta):
%             r = rc*(1 + sinh(beta*(xi - A))/sinh(beta*A)),  xi = (1+eta)/2
%   'malik' - r = a(1+eta)/(b-eta), половина точек в [0, rmap]
% =========================================================================
function [r, Dr, w] = radial_grid(N, g)
[De, eta] = cheb(N-1);               % eta убывает от 1 до -1
eta = flipud(eta); De = rot90(De, 2); % теперь eta возрастает
switch lower(g.mapping)
    case 'malik'
        if g.rmax <= 2*g.rmap, error('pse_jet:grid', 'Требуется rmax > 2*rmap'); end
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
        error('pse_jet:grid', 'Неизвестное отображение %s', g.mapping);
end
Dr = diag(1./drde)*De;
w  = clencurt(N-1);
w  = flipud(w(:)).*drde;             % веса для int f dr
r(1) = 0; r(end) = g.rmax;
end

function [D, x] = cheb(N)
% Trefethen, Spectral Methods in MATLAB
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x, 1, N+1);
dX = X - X';
D = (c*(1./c)')./(dX + eye(N+1));
D = D - diag(sum(D, 2));
end

function w = clencurt(N)
% Trefethen: веса квадратуры Кленшоу-Кертиса
theta = pi*(0:N)'/N;
w = zeros(1, N+1); ii = 2:N; v = ones(N-1, 1);
if mod(N, 2) == 0
    w(1) = 1/(N^2 - 1); w(N+1) = w(1);
    for k = 1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2 - 1); end
    v = v - cos(N*theta(ii))/(N^2 - 1);
else
    w(1) = 1/N^2; w(N+1) = w(1);
    for k = 1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2 - 1); end
end
w(ii) = 2*v/N;
end

% =========================================================================
%  Среднее течение
% =========================================================================
function mf = prepare_meanflow(mf)
mf.x = mf.x(:); mf.r = mf.r(:);
nx = numel(mf.x); nr = numel(mf.r);
chk = @(F, name) assert(isequal(size(F), [nx nr]), ...
    'pse_jet:mf', 'mf.%s должно иметь размер [numel(mf.x) x numel(mf.r)]', name);
chk(mf.U, 'U');
if ~isfield(mf, 'V') || isempty(mf.V), mf.V = zeros(nx, nr); end
chk(mf.V, 'V');
hasR = isfield(mf, 'rho') && ~isempty(mf.rho);
hasT = isfield(mf, 'T')   && ~isempty(mf.T);
if ~hasR && ~hasT, error('pse_jet:mf', 'Задайте mf.rho и/или mf.T'); end
if ~hasR, mf.rho = 1./mf.T; end
if ~hasT, mf.T = 1./mf.rho; end
chk(mf.rho, 'rho'); chk(mf.T, 'T');
% производные на исходной сетке (2-й порядок, неравномерная сетка)
[mf.Ur, mf.Ux] = grad2(mf.U,   mf.x, mf.r);
[mf.Vr, mf.Vx] = grad2(mf.V,   mf.x, mf.r);
[mf.rr, mf.rx] = grad2(mf.rho, mf.x, mf.r);
[mf.Tr, mf.Tx] = grad2(mf.T,   mf.x, mf.r);
end

function [Fr, Fx] = grad2(F, x, r)
Fx = zeros(size(F)); Fr = zeros(size(F));
if numel(x) > 1
    for j = 1:size(F,2), Fx(:,j) = d1(F(:,j), x); end
end
for i = 1:size(F,1), Fr(i,:) = d1(F(i,:).', r).'; end
end

function d = d1(f, x)
% производная 2-го порядка на неравномерной сетке
n = numel(x); d = zeros(n,1);
if n < 3, d(:) = (f(end)-f(1))/(x(end)-x(1)); return; end
h = diff(x(:));
for i = 2:n-1
    h1 = h(i-1); h2 = h(i);
    d(i) = (-h2/(h1*(h1+h2)))*f(i-1) + ((h2-h1)/(h1*h2))*f(i) + (h1/(h2*(h1+h2)))*f(i+1);
end
h1 = h(1); h2 = h(2);
d(1) = -(2*h1+h2)/(h1*(h1+h2))*f(1) + (h1+h2)/(h1*h2)*f(2) - h1/(h2*(h1+h2))*f(3);
h1 = h(n-1); h2 = h(n-2);
d(n) = (2*h1+h2)/(h1*(h1+h2))*f(n) - (h1+h2)/(h1*h2)*f(n-1) + h1/(h2*(h1+h2))*f(n-2);
end

function s = meanflow_at(mf, xq, r, method)
% Интерполяция среднего течения в сечение xq на радиальные точки r.
% За пределами исходной сетки по r значения "замораживаются", производные = 0.
xq = min(max(xq, mf.x(1)), mf.x(end));
rq = min(max(r, mf.r(1)), mf.r(end));
outside = r > mf.r(end) | r < mf.r(1);
names = {'U','V','rho','T','Ux','Ur','Vx','Vr','rx','rr','Tx','Tr'};
for k = 1:numel(names)
    F = mf.(names{k});
    if numel(mf.x) > 1
        row = interp1(mf.x, F, xq, method);
    else
        row = F;
    end
    v = interp1(mf.r, row(:), rq, method);
    if k > 4, v(outside) = 0; end
    s.(names{k}) = v(:);
end
s.src.mf = mf; s.src.x = xq;
end

function Fo = interp_r(r, F, ro)
Fo = interp1(r, real(F), ro, 'spline') + 1i*interp1(r, imag(F), ro, 'spline');
Fo(ro > r(end) | ro < r(1), :) = 0;
end

function opts = set_defaults(opts, mf, nSt)
def = struct('M', 0.9, 'gamma', 1.4, 'Nr', 120, 'rmax', 15, ...
    'mapping', 'sinh', 'rc', 0.5, 'beta', 5, 'rmap', 1.0, ...
    'x0', mf.x(1), 'x1', mf.x(end), 'dx', 0.1, 'dxAuto', true, 'dxSafety', 1.1, ...
    'dropDpdx', false, 'tol', 1e-8, 'maxIter', 30, 'normVars', 1:5, ...
    'alphaGuess', [], 'cphRange', [0.2 2], 'amp', 1, 'interp', 'pchip', ...
    'storeQ', false, 'rOut', [], 'verbose', true, 'lstTol', 1e-2, ...
    'minGrowth', 0.02);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || (isempty(opts.(f{i})) && ~isempty(def.(f{i})))
        opts.(f{i}) = def.(f{i});
    end
end
if isempty(opts.alphaGuess), opts.alphaGuess = nan(1, nSt); end
if isscalar(opts.alphaGuess), opts.alphaGuess = repmat(opts.alphaGuess, 1, nSt); end
end
