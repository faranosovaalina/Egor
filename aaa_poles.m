function p = aaa_poles(f)
%AAA_POLES  Полюса AAA-аппроксимации f (см. aaa_fit).

m = numel(f.wj);
B = eye(m + 1); B(1, 1) = 0;
E = [0, f.wj.'; ones(m, 1), diag(f.zj)];
p = eig(E, B);
p = p(isfinite(p));
end
