function [gs, gu, W] = pareto_span(a_s, h_u, nt, nph)
%PARETO_SPAN  Rate-sensing Pareto boundary of one beam toward two directions.
%   The optimal unit-norm beam lies in span{a_s, h_u}, so we search
%   w = cos(t) a_s/|a_s| + sin(t) e^{j phi} q   (q: h_u orthogonalised to a_s)
%   on a dense (t, phi) grid and keep, for each sensing gain, the largest
%   user gain. Returns gains |a_s' w|^2 and |h_u' w|^2 along the boundary.
if nargin < 3, nt = 401; end
if nargin < 4, nph = 72; end
u1 = a_s / norm(a_s);
q = h_u - u1 * (u1' * h_u);
if norm(q) < 1e-9, q = zeros(size(u1)); else, q = q / norm(q); end
t = linspace(0, pi/2, nt);  ph = linspace(0, 2*pi, nph + 1);  ph(end) = [];
[T, PH] = meshgrid(t, ph);
Wall = u1 * cos(T(:)).' + q * (sin(T(:)) .* exp(1j * PH(:))).';
Gs = abs(a_s' * Wall).^2;  Gu = abs(h_u' * Wall).^2;
edges = linspace(0, max(Gs), 300);
gs = zeros(1, numel(edges));  gu = gs;  W = zeros(numel(u1), numel(edges));
for k = 1:numel(edges)
    ok = find(Gs >= edges(k) - 1e-12);
    [gu(k), j] = max(Gu(ok));
    gs(k) = Gs(ok(j));  W(:, k) = Wall(:, ok(j));
end
end
