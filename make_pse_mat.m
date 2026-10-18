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
