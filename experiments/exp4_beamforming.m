function r = exp4_beamforming()
%EXP4_BEAMFORMING  One antenna array, two jobs: a beam for the user and one for the target.
%   32-element half-wavelength ULA. User at -35 deg, target at +25 deg.
%   The Pareto boundary is found exactly by searching the span of the two
%   steering vectors (PARETO_SPAN). Compared with three simple designs:
%   * power split: w ~ sqrt(rho) a_s + sqrt(1-rho) h_u (no phase alignment),
%   * weighted sum: principal eigenvector of rho a_s a_s' + (1-rho) h_u h_u',
%   * time sharing between the two pure beams.
c = palette();
Nt = 32;  th_u = -35;  th_s = 25;
h_u = ula(Nt, th_u);  a_s = ula(Nt, th_s);
snr_u0 = 1;                    % user SNR per unit beam gain: 0 dB -> 15 dB with full gain
rate = @(G) log2(1 + snr_u0 * G);
[cases, names] = deal({th_s - 60, th_s - 5}, {'60 deg apart', '5 deg apart'});
rho = linspace(0, 1, 101);
F = struct();
for k = 1:2
    hu = ula(Nt, cases{k});
    [gs, gu, Wf] = pareto_span(a_s, hu);
    Gps = zeros(2, numel(rho));  Gws = Gps;
    for i = 1:numel(rho)
        ws = sqrt(rho(i)) * a_s + sqrt(1 - rho(i)) * hu;  ws = ws / norm(ws);
        Gps(:, i) = [abs(a_s' * ws)^2; abs(hu' * ws)^2];
        A = rho(i) * (a_s * a_s') + (1 - rho(i)) * (hu * hu');
        [V, D] = eig((A + A') / 2);  [~, j] = max(real(diag(D)));  w = V(:, j);
        Gws(:, i) = [abs(a_s' * w)^2; abs(hu' * w)^2];
    end
    F(k).gs = gs;  F(k).gu = gu;  F(k).W = Wf;  F(k).ps = Gps;  F(k).ws = Gws;
end
% separation sweep: sensing gain left when the user keeps 90 % of its best rate
seps = [1 2 3 4 6 8 10 15 20 40 60];
keep = zeros(size(seps));
for j = 1:numel(seps)
    hu = ula(Nt, th_s - seps(j));
    [gs, gu] = pareto_span(a_s, hu);
    ok = rate(gu) >= 0.9 * rate(Nt);
    keep(j) = max(gs(ok)) / Nt;
end
h = new_figure(1750, 620);
ax = subplot(1, 3, 1);  hold on;
col = {c.teal, c.violet};
for k = 1:2
    hh(1) = plot(10*log10(F(k).gs/Nt), rate(F(k).gu), '-', 'color', col{k}, 'linewidth', 3.2);
    hh(2) = plot(10*log10(F(k).ps(1, :)/Nt), rate(F(k).ps(2, :)), '--', 'color', c.gold, 'linewidth', 1.8);
    hh(3) = plot(10*log10(F(k).ws(1, :)/Nt), rate(F(k).ws(2, :)), 'o', 'color', c.coral, 'markersize', 4);
end
hh(4) = plot(10*log10([F(1).gs(1) F(1).gs(end)]/Nt), rate([F(1).gu(1) F(1).gu(end)]), ':', 'color', c.muted, 'linewidth', 2);
text(-26, 4.75, names{1}, 'color', c.teal, 'fontweight', 'bold');
text(-7.5, 3.0, names{2}, 'color', c.violet, 'fontweight', 'bold');
xlim([-30 1]);  ylim([0 5.5]);
xlabel('sensing beam gain at the target (dB re. full beam)');
ylabel('user rate (bit/s/Hz)');
title('Rate vs sensing: the Pareto boundary');
style_axes(ax, c);
lg = legend(hh, {'Pareto boundary (span search)', 'power split', 'weighted-sum eigenbeam', 'time sharing (60 deg)'}, 'location', 'southwest');
set(lg, 'color', c.bg, 'textcolor', c.fg, 'edgecolor', c.grid, 'fontsize', 10);
ax = subplot(1, 3, 2);  hold on;
th = linspace(-90, 90, 721);  A = ula(Nt, th);
[~, i1] = min(abs(10*log10(F(1).gs/Nt) - (-3)));
picks = {F(1).W(:, 1), F(1).W(:, i1), F(1).W(:, end)};
cols = {c.coral, c.gold, c.teal};
for k = 1:3
    g = 10*log10(abs(A' * picks{k}).^2);
    plot(th, max(g, -15), 'color', cols{k}, 'linewidth', 2);
end
plot([th_u th_u], [-15 17], ':', 'color', c.muted);  plot([th_s th_s], [-15 17], ':', 'color', c.muted);
text(th_u - 13, 17, 'user', 'color', c.fg);  text(th_s + 2, 17, 'target', 'color', c.fg);
xlim([-90 90]);  ylim([-15 18.5]);
xlabel('angle (deg)');  ylabel('array gain (dBi)');
title('Beam patterns along the boundary');
style_axes(ax, c);
dark_legend(ax, {'all to the user', 'target at -3 dB', 'all to the target'}, 'south');
ax = subplot(1, 3, 3);
bar(ax, 10*log10(keep), 'facecolor', c.violet);
set(ax, 'xticklabel', arrayfun(@(s) sprintf('%d', s), seps, 'uniformoutput', false));
xlabel('angle between user and target (deg)');
ylabel('sensing gain left (dB re. full beam)');
title('Sensing left when the user keeps 90 % of its rate');
ylim([min(10*log10(keep)) - 1.5 0.5]);
style_axes(ax, c);
for j = 1:numel(seps)
    text(j - 0.4, 10*log10(keep(j)) - 0.6, sprintf('%.1f', 10*log10(keep(j))), 'color', c.fg, 'fontsize', 9);
end
save_figure(h, 'fig4_beamforming_pareto');
write_csv('exp4_pareto_60deg', 'target_gain,user_gain,rate_bpsHz', [F(1).gs.' F(1).gu.' rate(F(1).gu).']);
write_csv('exp4_pareto_5deg', 'target_gain,user_gain,rate_bpsHz', [F(2).gs.' F(2).gu.' rate(F(2).gu).']);
write_csv('exp4_separation', 'separation_deg,sensing_gain_kept_dB', [seps.' 10*log10(keep).']);
% headline numbers: user rate when the target gets half of its full beam (-3 dB)
for k = 1:2
    [~, i] = min(abs(10*log10(F(k).gs/Nt) + 3));
    r.front3(k) = rate(F(k).gu(i));
    ok = F(k).ps(1, :) >= F(k).gs(i);  r.split3(k) = max(rate(F(k).ps(2, ok)));
    ok = F(k).ws(1, :) >= F(k).gs(i);  r.ws3(k) = max([0 rate(F(k).ws(2, ok))]);
end
r.ends = [rate(Nt) rate(F(1).gu(end))];
r.keep_db = 10*log10(keep);  r.seps = seps;
end
