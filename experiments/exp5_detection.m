function r = exp5_detection()
%EXP5_DETECTION  Probability of detection with CFAR, and how far each object is seen.
%   Single target on a bin centre, division receiver, QPSK, Hann windows.
%   Theory (known noise, non-fluctuating target):
%   Pd = Q1( sqrt(2 SNR_int), sqrt(-2 ln Pfa) ).
p = isac_params();  c = palette();
randn('seed', 5);  rand('seed', 5);
pfa = 1e-4;  guard = [2 2];  train = [8 6];
win_loss = 10*log10((0.5^2 / 0.375)^2);           % Hann coherent loss, both dims (-3.5 dB)
gain = 10*log10(p.N * p.M) + win_loss;
snr_int = 6:1:20;                                 % integrated SNR after windows (dB)
trials = 120;
Rt = 40 * p.dR;  vt = 10 * p.dv;
pd_ca = zeros(size(snr_int));  pd_so = pd_ca;
for i = 1:numel(snr_int)
    a = 10^((snr_int(i) - gain) / 20);
    hit_ca = 0;  hit_so = 0;
    for t = 1:trials
        tg = struct('R', Rt, 'v', vt, 'a', a * exp(2j*pi*rand));
        X = qam_mod(4, [p.N p.M]);
        Y = X .* radar_channel(p, tg) + (randn(p.N, p.M) + 1j*randn(p.N, p.M)) / sqrt(2);
        [P, rax, vax] = rd_map(p, Y, X, 'div', 1, true);
        [~, ir] = min(abs(rax - Rt));  [~, iv] = min(abs(vax - vt));
        rows = mod(ir - 12 + (0:22), p.N) + 1;  cols = mod(iv - 10 + (0:18), p.M) + 1;   % target at (12, 10)
        sub = P(rows, cols);
        d1 = cfar2d(sub, guard, train, pfa, 'ca');  d2 = cfar2d(sub, guard, train, pfa, 'so');
        hit_ca = hit_ca + d1(12, 10);  hit_so = hit_so + d2(12, 10);
    end
    pd_ca(i) = hit_ca / trials;  pd_so(i) = hit_so / trials;
end
pd_th = marcum_q1(sqrt(2 * 10.^(snr_int/10)), sqrt(-2 * log(pfa)));
% measured false-alarm rate on noise-only maps
fa = [0 0];  cells = 0;
for t = 1:12
    X = qam_mod(4, [p.N p.M]);
    W = (randn(p.N, p.M) + 1j*randn(p.N, p.M)) / sqrt(2);
    P = rd_map(p, W, X, 'div', 1, true);
    fa(1) = fa(1) + sum(sum(cfar2d(P, guard, train, pfa, 'ca')));
    fa(2) = fa(2) + sum(sum(cfar2d(P, guard, train, pfa, 'so')));
    cells = cells + numel(P);
end
pfa_meas = fa / cells;
% detection range for Pd = 0.9 at Pfa = 1e-6 (theory) for typical road users
need = 0;
for s = 0:0.01:25
    if marcum_q1(sqrt(2*10^(s/10)), sqrt(-2*log(1e-6))) >= 0.9, need = s; break; end
end
names = {'pedestrian', 'cyclist', 'car', 'truck'};  rcs = [0.5 1 10 100];
Pt = [20 30 40];
Rdet = zeros(numel(rcs), numel(Pt));
for k = 1:numel(rcs)
    for j = 1:numel(Pt)
        q = p;  q.Pt_dBm = Pt(j);
        snr1 = echo_snr(q, 1, rcs(k));                  % at 1 m
        Rdet(k, j) = (snr1 / 10^((need - gain)/10))^(1/4);
    end
end
h = new_figure(1600, 620);
ax = subplot(1, 2, 1);
plot(snr_int, pd_th, '-', 'color', c.teal, 'linewidth', 3);  hold on;
plot(snr_int, pd_ca, 'o-', 'color', c.coral, 'markerfacecolor', c.coral, 'linewidth', 1.6);
plot(snr_int, pd_so, 's-', 'color', c.gold, 'markerfacecolor', c.gold, 'linewidth', 1.6);
xlabel('integrated SNR in the target cell (dB)');  ylabel('probability of detection');
title(sprintf('P_{fa} target 10^{-4}: measured CA %.1e, SO %.1e', pfa_meas));
ylim([0 1.02]);
style_axes(ax, c);
dark_legend(ax, {'theory, known noise', sprintf('CA-CFAR, %d training cells', prod(2*(guard+train)+1) - prod(2*guard+1)), 'SO-CFAR'}, 'southeast');
ax = subplot(1, 2, 2);
b = barh(ax, Rdet);
cols = {c.sky, c.teal, c.coral};
for j = 1:numel(Pt), set(b(j), 'facecolor', cols{j}); end
hold on;
plot([p.Rcp p.Rcp], [0.4 numel(rcs) + 0.6], '--', 'color', c.gold, 'linewidth', 2);
text(p.Rcp * 1.03, 0.55, sprintf('cyclic prefix limit %.0f m', p.Rcp), 'color', c.gold, 'fontsize', 10);
set(ax, 'yticklabel', names, 'xscale', 'log');
xlim([30 1e4]);
xlabel('range for P_d = 0.9, P_{fa} = 10^{-6} (m)');
title(sprintf('How far each road user is seen (one %.1f ms frame)', p.M * p.Ts * 1e3));
style_axes(ax, c);
dark_legend(ax, arrayfun(@(x) sprintf('%d dBm', x), Pt, 'uniformoutput', false), 'southeast');
save_figure(h, 'fig5_detection');
write_csv('exp5_pd', 'snr_int_dB,pd_theory,pd_ca_cfar,pd_so_cfar', [snr_int.' pd_th.' pd_ca.' pd_so.']);
write_csv('exp5_range', 'rcs_m2,range_20dBm,range_30dBm,range_40dBm', [rcs.' Rdet]);
r.pfa_meas = pfa_meas;  r.need = need;  r.Rdet = Rdet;
r.loss = [xcross(snr_int, pd_th, 0.9) xcross(snr_int, pd_ca, 0.9) xcross(snr_int, pd_so, 0.9)];
end

function x = xcross(xs, ys, y0)
% first crossing of y0 by linear interpolation
k = find(ys >= y0, 1);
if isempty(k) || k == 1, x = NaN; return; end
x = xs(k-1) + (y0 - ys(k-1)) * (xs(k) - xs(k-1)) / (ys(k) - ys(k-1));
end
