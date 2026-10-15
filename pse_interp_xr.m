function v = pse_interp_xr(x, r, F, xq, rq)
%PSE_INTERP_XR  Интерполяция комплексного поля F(r,x) с прямоугольной
%сетки PSE в произвольные точки (xq, rq).
%   F: [numel(r) x numel(x)];  xq, rq - массивы одного размера.
% Сплайн по x, затем кубический Лагранж по r
% (interp2 'spline' в разрозненных точках в Octave не работает).
Fx = interp1(x(:), real(F.'), xq(:), 'spline') + ...
     1i*interp1(x(:), imag(F.'), xq(:), 'spline');      % [Nq x Nr]
if numel(xq) == 1, Fx = Fx(:).'; end
r  = r(:); n = numel(r); rq = rq(:);
i0 = zeros(size(rq));
for q = 1:numel(rq)
    i0(q) = find(r <= rq(q), 1, 'last');
end
i0 = min(max(i0 - 1, 1), n - 3);            % опорные точки i0..i0+3
v  = zeros(size(rq));
for q = 1:numel(rq)
    idx = i0(q) + (0:3);
    rr  = r(idx);
    L   = ones(1, 4);
    for a = 1:4
        for b = [1:a-1, a+1:4]
            L(a) = L(a)*(rq(q) - rr(b))/(rr(a) - rr(b));
        end
    end
    v(q) = Fx(q, idx)*L.';
end
v = reshape(v, size(xq));
end
