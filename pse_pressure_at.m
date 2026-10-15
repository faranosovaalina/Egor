function p = pse_pressure_at(src, iSt, pts)
%PSE_PRESSURE_AT  Давление PSE (прямой расчёт) в заданных точках (x, r) -
%для сравнения с результатом kirchhoff_pse в ближнем поле.
%
%   p = pse_pressure_at(src, iSt, pts)
%
%   src - mat-файл или структура (формат как в kirchhoff_pse)
%   iSt - номер числа Струхаля в файле
%   pts - [Npts x 2] точки (x, r) внутри PSE-сетки
%   p   - [Npts x 1] комплексное давление (плоскость phi = 0)

if ischar(src), S = load(src); else S = src; end
xa = S.x_all;     if iscell(xa), xa = xa{iSt}; end
P  = S.P_all_St;  if iscell(P),  P  = P{iSt};  end
x  = xa(:).';
r  = S.r_grid(:);
if isfield(S, 'alpha_all') && ~isempty(S.alpha_all)
    al = S.alpha_all; if iscell(al), al = al{iSt}; end
    al = al(:).';
else
    al = zeros(size(x));
end
Phi = cumtrapz(x, al);
H   = P.*repmat(exp(-1i*Phi), numel(r), 1);
Hq  = pse_interp_xr(x, r, H, pts(:,1), pts(:,2));
Pq  = interp1(x, real(Phi), pts(:,1), 'spline') + 1i*interp1(x, imag(Phi), pts(:,1), 'spline');
p   = Hq(:).*exp(1i*Pq(:));
end
