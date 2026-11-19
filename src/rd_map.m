function [P, r_ax, v_ax] = rd_map(p, Y, X, method, zp, win)
%RD_MAP  Range-Doppler periodogram of received symbols Y given sent symbols X.
%   method 'div'  : Z = Y ./ X  (removes the data, Sturm & Wiesbeck)
%          'corr' : Z = Y .* conj(X) (matched filter, keeps the data in the sidelobes)
%   zp: zero padding factor, win: true for Hann windows.
%   P is |.|^2 normalised so that a unit-amplitude echo in white noise of
%   variance s2 gives a peak of N*M (coherent gain) over a floor of s2.
if nargin < 4, method = 'div'; end
if nargin < 5, zp = 1; end
if nargin < 6, win = true; end
if strcmp(method, 'div'), Z = Y ./ X; else, Z = Y .* conj(X); end
if win
    wn = hann_win(p.N);  wm = hann_win(p.M).';
    Z = Z .* (wn * wm) / sqrt(mean(wn.^2) * mean(wm.^2));
end
Nr = zp * p.N;  Nd = zp * p.M;
r = ifft(Z, Nr, 1) * Nr;                % range profile per symbol (sum over n)
P = abs(fftshift(fft(r, Nd, 2), 2)).^2 / (p.N * p.M);
r_ax = (0:Nr-1) * p.c / (2 * p.df * Nr);
fD = ((0:Nd-1) - floor(Nd/2)) / (Nd * p.Ts);
v_ax = -fD * p.lambda / 2;
end
