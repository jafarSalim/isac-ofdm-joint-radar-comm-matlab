function H = radar_channel(p, tg)
%RADAR_CHANNEL  Frequency-domain echo channel H(n, m) of point targets.
%   tg: struct array with fields R (m), v (m/s, positive = approaching? no:
%   positive = moving away) and a (complex amplitude, linear).
%   H(n,m) = sum_k a_k exp(-j 2 pi n df tau_k) exp(j 2 pi m Ts fD_k),
%   tau = 2R/c, fD = -2v/lambda. Valid while tau < Tcp and fD << df.
n = (0:p.N-1).';  m = 0:p.M-1;
H = zeros(p.N, p.M);
for k = 1:numel(tg)
    tau = 2 * tg(k).R / p.c;
    fD  = -2 * tg(k).v / p.lambda;
    H = H + tg(k).a * exp(-2j*pi*n*p.df*tau) * exp(2j*pi*m*p.Ts*fD);
end
end
