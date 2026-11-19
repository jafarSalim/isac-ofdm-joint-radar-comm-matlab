function r = exp3_modulation()
%EXP3_MODULATION  What the user data does to the radar, and what it gives the user.
%   (a) data-induced floor of the correlation receiver vs the noise-enhancement
%       loss of the division receiver, for QPSK, 16-QAM, 64-QAM;
%   (b) bit error rate of the same frame decoded by the user (AWGN, perfect CSI).
p = isac_params();  c = palette();
randn('seed', 3);  rand('seed', 3);
mods = [4 16 64];  names = {'QPSK', '16-QAM', '64-QAM'};
tg = struct('R', 30, 'v', 5, 'a', 10^(20/20));      % strong echo, 20 dB per element
H = radar_channel(p, tg);
tab = zeros(3, 6);
for k = 1:3
    X = qam_mod(mods(k), [p.N p.M]);
    W = (randn(p.N, p.M) + 1j*randn(p.N, p.M)) / sqrt(2);
    for meth = 1:2
        if meth == 1, mt = 'div'; else, mt = 'corr'; end
        [P, rax] = rd_map(p, X .* H + W, X, mt, 1, true);
        far = rax > 200 & rax < 1000;                 % far from the target
        floor_db = 10*log10(mean(mean(P(far, :))));
        tab(k, meth) = floor_db;
        tab(k, 2 + meth) = 10*log10(max(P(:))) - floor_db;   % peak to floor
    end
    tab(k, 5) = 10*log10(mean(1 ./ abs(X(:)).^2));            % theory: division loss
    tab(k, 6) = 10*log10(var(abs(X(:)).^2) * abs(tg.a)^2 + 1); % theory: correlation floor
end
% communication: BER of the same OFDM symbols over AWGN
es = 0:1:30;  ber = zeros(3, numel(es));
for k = 1:3
    for i = 1:numel(es)
        [X, b] = qam_mod(mods(k), [p.N 8]);
        Y = X + 10^(-es(i)/20) * (randn(size(X)) + 1j*randn(size(X))) / sqrt(2);
        ber(k, i) = mean(qam_demod(mods(k), Y) ~= b);
    end
end
h = new_figure(1700, 600);
ax = subplot(1, 3, 1);
b = bar(ax, tab(:, 1:2));  set(b(1), 'facecolor', c.teal);  set(b(2), 'facecolor', c.coral);
set(ax, 'xticklabel', names);
ylabel('noise floor of the map (dB, noise = 0)');
title('Where the floor ends up');
style_axes(ax, c);  dark_legend(ax, {'division', 'correlation'}, 'northwest');
for k = 1:3
    text(k - 0.3, tab(k, 1) + 1.2, sprintf('%.1f', tab(k, 1)), 'color', c.teal, 'fontsize', 10, 'fontweight', 'bold');
    text(k + 0.06, tab(k, 2) + 1.2, sprintf('%.1f', tab(k, 2)), 'color', c.coral, 'fontsize', 10, 'fontweight', 'bold');
end
ylim([0 max(tab(:, 2)) + 6]);
ax = subplot(1, 3, 2);
b = bar(ax, tab(:, 3:4));  set(b(1), 'facecolor', c.teal);  set(b(2), 'facecolor', c.coral);
set(ax, 'xticklabel', names);
ylabel('peak / floor (dB)');
title(sprintf('Dynamic range for a %.0f dB echo', 20 + 10*log10(p.N*p.M)));
style_axes(ax, c);  dark_legend(ax, {'division', 'correlation'}, 'southwest');
ylim([0 85]);
ax = subplot(1, 3, 3);
for k = 1:3
    semilogy(es, max(ber(k, :), 1e-7), '-o', 'color', c.mods{k}, 'markerfacecolor', c.mods{k}, ...
        'markersize', 4, 'linewidth', 2);  hold on;
end
ylim([1e-6 1]);  xlabel('user SNR (dB)');  ylabel('bit error rate');
title('What the user gets: 2, 4 or 6 bits per element');
style_axes(ax, c);  dark_legend(ax, names, 'southwest');
save_figure(h, 'fig3_modulation_tradeoff');
write_csv('exp3_radar', 'M,floor_div_dB,floor_corr_dB,peak_to_floor_div_dB,peak_to_floor_corr_dB,theory_div_loss_dB,theory_corr_floor_dB', [mods.' tab]);
write_csv('exp3_ber', 'Es_N0_dB,ber_qpsk,ber_16qam,ber_64qam', [es.' ber.']);
r.tab = tab;
for k = 1:3
    r.snr_1e3(k) = es(find(ber(k, :) < 1e-3, 1));
end
r.rate_gbps = log2(mods) * p.N * p.M / (p.M * p.Ts) / 1e9;
end
