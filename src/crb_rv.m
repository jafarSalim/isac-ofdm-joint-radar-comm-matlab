function [sR, sv] = crb_rv(p, snr)
%CRB_RV  Cramer-Rao bounds (standard deviations) of range and velocity.
%   One target, unknown complex amplitude, known data symbols of unit
%   average energy, per-element SNR snr (linear), N subcarriers, M symbols:
%   var(tau) >= 1 / (2 snr (2 pi df)^2 M N (N^2 - 1)/12)
%   var(fD)  >= 1 / (2 snr (2 pi Ts)^2 N M (M^2 - 1)/12)
Jt = 2 * snr * (2*pi*p.df)^2 * p.M * p.N * (p.N^2 - 1) / 12;
Jf = 2 * snr * (2*pi*p.Ts)^2 * p.N * p.M * (p.M^2 - 1) / 12;
sR = p.c / 2 ./ sqrt(Jt);
sv = p.lambda / 2 ./ sqrt(Jf);
end
