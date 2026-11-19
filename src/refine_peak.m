function [R, v] = refine_peak(p, Z, R0, v0)
%REFINE_PEAK  Maximum-likelihood refinement of one target around (R0, v0).
%   Maximises |sum Z(n,m) exp(+j2pi n df tau) exp(-j2pi m Ts fD)| by
%   alternating parabolic steps on a shrinking grid (no toolbox).
n = (0:p.N-1).';  m = (0:p.M-1);
R = R0;  v = v0;
sR = p.dR / 4;  sv = p.dv / 4;
for it = 1:7
    R = step1(@(x) cost(x, v), R, sR);
    v = step1(@(x) cost(R, x), v, sv);
    sR = sR / 3;  sv = sv / 3;
end
    function f = cost(Rx, vx)
        tau = 2 * Rx / p.c;  fD = -2 * vx / p.lambda;
        f = abs(exp(2j*pi*n*p.df*tau).' * Z * exp(-2j*pi*m*p.Ts*fD).');
    end
end

function x = step1(f, x0, h)
y = [f(x0 - h), f(x0), f(x0 + h)];
d = y(1) - 2*y(2) + y(3);
if d < 0
    x = x0 + 0.5 * h * (y(1) - y(3)) / d;
    x = min(max(x, x0 - h), x0 + h);
else
    [~, i] = max(y);  x = x0 + (i - 2) * h;
end
end
