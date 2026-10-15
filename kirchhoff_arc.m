function obs = kirchhoff_arc(Rarc, thetaDeg, xc)
%KIRCHHOFF_ARC  Точки на полярной дуге в плоскости (x, r) для kirchhoff_pse.
%
%   obs = kirchhoff_arc(Rarc, thetaDeg)
%   obs = kirchhoff_arc(Rarc, thetaDeg, xc)
%
%   Rarc     - радиус дуги (в D; например 200*R = 100 при R = 0.5)
%   thetaDeg - полярные углы от оси струи (вниз по потоку), градусы
%   xc       - x центра дуги (по умолчанию 0 - срез сопла)
%   obs      - [numel(thetaDeg) x 2] = [x, r]

if nargin < 3, xc = 0; end
th  = thetaDeg(:)*pi/180;
obs = [xc + Rarc*cos(th), Rarc*sin(th)];
end
