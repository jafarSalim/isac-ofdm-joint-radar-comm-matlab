function r = exp2_crb()
%EXP2_CRB  Range and velocity accuracy against the Cramer-Rao bound.
%   One target (R = 50.3 m, v = 7.7 m/s), QPSK data, division receiver,
%   coarse peak of an unwindowed 2-D FFT, then maximum-likelihood refinement.
p = isac_params();  c = palette();
randn('seed', 2);  rand('seed', 2);
snr_db = -50:2.5:-10;  trials = 60;
R0 = 50.3;  v0 = 7.7;
eR = zeros(numel(snr_db), trials);  ev = eR;
n = (0:p.N-1).';  m = 0:p.M-1;
for i = 1:numel(snr_db)
    a = 10^(snr_db(i)/20);
    for t = 1:trials
        tg = struct('R', R0, 'v', v0, 'a', a * exp(2j*pi*rand));
        X = qam_mod(4, [p.N p.M]);
        Y = X .* radar_channel(p, tg) + (randn(p.N, p.M) + 1j*randn(p.N, p.M)) / sqrt(2);
        Z = Y ./ X;
        [P, rax, vax] = rd_map(p, Y, X, 'div', 2, false);
        P(rax > 150, :) = 0;                       % search window 0-150 m
        [~, k] = max(P(:));  [ir, iv] = ind2sub(size(P), k);
        [Rh, vh] = refine_peak(p, Z, rax(ir), vax(iv));
        eR(i, t) = Rh - R0;  ev(i, t) = vh - v0;
    end
end
rmseR = sqrt(mean(eR.^2, 2));  rmsev = sqrt(mean(ev.^2, 2));
[cR, cv] = crb_rv(p, 10.^(snr_db/10));
outR = mean(abs(eR) > 3 * p.dR, 2);  outv = mean(abs(ev) > 3 * p.dv, 2);
gain = 10*log10(p.N * p.M);
h = new_figure(1600, 620);
lab = {'range (m)', 'velocity (m/s)'};
dat = {rmseR, rmsev; cR, cv; outR, outv};
res = [p.dR p.dv];
for k = 1:2
    ax = subplot(1, 2, k);
    semilogy(snr_db, dat{2, k}, '-', 'color', c.teal, 'linewidth', c.lw);  hold on;
    semilogy(snr_db, dat{1, k}, 'o', 'color', c.coral, 'markerfacecolor', c.coral, 'markersize', 7);
    semilogy(snr_db, res(k) / sqrt(12) * ones(size(snr_db)), '--', 'color', c.gold, 'linewidth', 1.6);
    xlabel(sprintf('SNR per resource element (dB);  integrated = +%.1f dB', gain));
    ylabel(['RMSE of ' lab{k}]);
    un = {'m', 'm/s'};
    title(sprintf('%s: resolution %.2f %s, RMSE at -20 dB %.3g %s', strtok(lab{k}), res(k), un{k}, dat{1, k}(snr_db == -20), un{k}));
    style_axes(ax, c);
    dark_legend(ax, {'Cramer-Rao bound', 'FFT peak + ML refinement', 'FFT bin only (resolution/sqrt 12)'}, 'northeast');
end
save_figure(h, 'fig2_accuracy_vs_crb');
write_csv('exp2_crb', 'snr_dB,rmse_range_m,crb_range_m,rmse_velocity_mps,crb_velocity_mps,outlier_range,outlier_velocity', ...
    [snr_db.' rmseR cR.' rmsev cv.' outR outv]);
r.snr_db = snr_db;  r.ratioR = rmseR ./ cR.';  r.ratiov = rmsev ./ cv.';
hi = snr_db >= -35;
r.eff = [mean(r.ratioR(hi)) mean(r.ratiov(hi))];
r.threshold = snr_db(find(outR < 0.05, 1));
r.at20 = [rmseR(snr_db == -20) rmsev(snr_db == -20)];
end
