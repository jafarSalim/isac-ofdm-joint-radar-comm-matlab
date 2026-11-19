function snr = echo_snr(p, R, rcs)
%ECHO_SNR  Per-resource-element SNR (linear) of a radar echo from the radar equation.
%   Pr = Pt G^2 lambda^2 rcs / ((4 pi)^3 R^4); noise kT F B over the band.
Pt = 10^((p.Pt_dBm - 30) / 10);  G = 10^(p.G_dBi / 10);
Pr = Pt * G^2 * p.lambda^2 * rcs ./ ((4*pi)^3 * R.^4);
N0B = 10^((p.kT_dBm + p.NF_dB - 30) / 10) * p.B;
snr = Pr ./ N0B;
end
