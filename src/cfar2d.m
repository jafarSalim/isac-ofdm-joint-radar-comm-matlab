function [det, thr] = cfar2d(P, guard, train, pfa, mode)
%CFAR2D  CFAR detector on a 2-D power map (square-law, exponential noise).
%   guard, train: half sizes [rows cols] (rows = range, cols = Doppler).
%   mode 'ca' (default): cell averaging over the whole training ring,
%        threshold = alpha * mean, alpha = Ntr (pfa^(-1/Ntr) - 1).
%   mode 'so': smallest-of the two halves in range (nearer / farther cells).
%        A strong neighbour on one side no longer masks the cell under test.
%        alpha uses Ntr/2 cells (the half that is averaged).
if nargin < 5, mode = 'ca'; end
g = guard;  t = guard + train;
K = ones(2*t + 1);
K(t(1)+1-g(1):t(1)+1+g(1), t(2)+1-g(2):t(2)+1+g(2)) = 0;
Pp = [P(end-t(1)+1:end, :); P; P(1:t(1), :)];          % range and Doppler wrap
Pp = [Pp(:, end-t(2)+1:end), Pp, Pp(:, 1:t(2))];
if strcmp(mode, 'ca')
    Ntr = sum(K(:));
    S = conv2(Pp, K, 'valid') / Ntr;
else
    Ka = K;  Ka(t(1)+1:end, :) = 0;                       % rows before the cell
    Kb = K;  Kb(1:t(1)+1, :) = 0;                         % rows after the cell
    Ntr = sum(Ka(:));
    % conv2 flips the kernel; flip back so 'before' really means before
    Sa = conv2(Pp, rot90(Ka, 2), 'valid') / Ntr;
    Sb = conv2(Pp, rot90(Kb, 2), 'valid') / Ntr;
    S = min(Sa, Sb);
end
alpha = Ntr * (pfa^(-1/Ntr) - 1);
thr = alpha * S;
det = P > thr;
end
