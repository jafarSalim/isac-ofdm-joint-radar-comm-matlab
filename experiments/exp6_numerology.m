function r = exp6_numerology()
%EXP6_NUMEROLOGY  The cyclic prefix sets the sensing range, the numerology sets the rest.
%   (a) Time-domain simulation (numerology 3): a target whose echo arrives
%       after the cyclic prefix causes inter-symbol and inter-carrier
%       interference. We measure the loss of the range-Doppler peak.
%   (b) Range, resolution and velocity limits of the four 5G numerologies.
c = palette();
randn('seed', 6);  rand('seed', 6);
p = isac_params(3);  p.M = 64;                       % shorter frame for the time-domain run
Rs = [20 40 60 80 88 95 110 130 160 200 250 300];
loss = zeros(size(Rs));  sinr = loss;
fs = p.N * p.df;
for i = 1:numel(Rs)
    X = qam_mod(4, [p.N p.M]);
    x = ifft(X, [], 1) * sqrt(p.N);                  % unit-power time samples
    x = [x(end-p.Ncp+1:end, :); x];                  % add CP
    s = x(:);
    d = round(2 * Rs(i) / p.c * fs);                 % echo delay in samples
    fD = 2 / (p.M * p.Ts);                           % on a Doppler bin (about 19 m/s)
    e = [zeros(d, 1); s(1:end-d)] .* exp(2j*pi*fD*(0:numel(s)-1).'/fs);
    E = reshape(e, p.N + p.Ncp, p.M);
    Y = fft(E(p.Ncp+1:end, :), [], 1) / sqrt(p.N);   % remove CP, back to subcarriers
    P = rd_map(p, Y, X, 'div', 1, false);
    [mx, k] = max(P(:));
    ideal = p.N * p.M;                               % noise-free peak without ISI
    loss(i) = 10*log10(ideal / mx);
    Pm = P;  Pm(k) = 0;
    sinr(i) = 10*log10(mx / mean(Pm(:)));            % peak over the interference floor
end
mus = 0:3;
tab = zeros(numel(mus), 6);
for j = 1:numel(mus)
    q = isac_params(mus(j));
    tab(j, :) = [q.df/1e3, q.B/1e6, q.dR, q.Rcp, q.dv, q.vmax];
end
h = new_figure(1600, 620);
ax = subplot(1, 2, 1);
plot(Rs, loss, 'o-', 'color', c.coral, 'markerfacecolor', c.coral, 'linewidth', 2.2);  hold on;
plot(Rs, -sinr + max(sinr), 's-', 'color', c.teal, 'markerfacecolor', c.teal, 'linewidth', 2.2);
yl = [0 max([loss -sinr+max(sinr)]) + 3];
plot([p.Rcp p.Rcp], yl, '--', 'color', c.gold, 'linewidth', 2);
text(p.Rcp - 4, yl(2) * 0.45, sprintf('echo leaves the CP\nat %.1f m', p.Rcp), 'color', c.gold, 'horizontalalignment', 'right', 'fontweight', 'bold');
ylim(yl);
xlabel('target range (m)');  ylabel('dB');
title('Numerology 3 (120 kHz): time-domain simulation');
style_axes(ax, c);
dark_legend(ax, {'peak loss', 'drop of peak-to-interference ratio'}, 'northwest');
ax = subplot(1, 2, 2);
semilogy(mus, tab(:, 4), 'o-', 'color', c.gold, 'markerfacecolor', c.gold, 'linewidth', 2.4);  hold on;
semilogy(mus, tab(:, 3), 's-', 'color', c.teal, 'markerfacecolor', c.teal, 'linewidth', 2.4);
semilogy(mus, tab(:, 6), '^-', 'color', c.coral, 'markerfacecolor', c.coral, 'linewidth', 2.4);
semilogy(mus, tab(:, 5), 'd-', 'color', c.violet, 'markerfacecolor', c.violet, 'linewidth', 2.4);
set(ax, 'xtick', mus, 'xticklabel', arrayfun(@(j) sprintf('\\mu=%d\n%g kHz', mus(j), tab(j,1)), 1:4, 'uniformoutput', false));
xlim([-0.3 3.3]);
ylabel('m  or  m/s');
title('Same 1024 subcarriers, four numerologies');
style_axes(ax, c);
dark_legend(ax, {'max range (CP), m', 'range resolution, m', 'max velocity, m/s', 'velocity resolution, m/s'}, 'east');
save_figure(h, 'fig6_cp_and_numerology');
write_csv('exp6_cp', 'range_m,peak_loss_dB,peak_to_interference_dB', [Rs.' loss.' sinr.']);
write_csv('exp6_numerology', 'scs_kHz,bandwidth_MHz,range_res_m,cp_range_m,velocity_res_mps,max_velocity_mps', tab);
r.loss = loss;  r.Rs = Rs;  r.tab = tab;  r.sinr = sinr;
end
