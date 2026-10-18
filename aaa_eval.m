function v = aaa_eval(f, z)
%AAA_EVAL  Значения AAA-аппроксимации f (см. aaa_fit) в точках z (комплексные
%допускаются); форма результата как у z.

sz = size(z);
z  = z(:);
Cz = 1./bsxfun(@minus, z, f.zj.');
v  = (Cz*(f.wj.*f.fj))./(Cz*f.wj);
bad = find(~isfinite(v));                 % точки, совпавшие с узлами
for q = bad.'
    [~, k] = min(abs(z(q) - f.zj));
    v(q) = f.fj(k);
end
v = reshape(v, sz);
end
