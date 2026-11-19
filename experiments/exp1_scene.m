function r = exp1_scene()
%EXP1_SCENE  Range-Doppler maps of a street scene from one OFDM data frame.
%   The frame carries real 16-QAM / 64-QAM user data. Two receivers:
%   element-wise division (data removed) and correlation (matched filter).
p = isac_params();  c = palette();
randn('seed', 1);  rand('seed', 1);
tg = scene(p, 0);
H = radar_channel(p, tg);
W = (randn(p.N, p.M) + 1j * randn(p.N, p.M)) / sqrt(2);
cases = {'div', 16; 'corr', 64};
tit = {'Division receiver, 16-QAM data', 'Correlation receiver, 64-QAM data'};
h = new_figure(1600, 980);
r.rows = [];
for k = 1:2
    X = qam_mod(cases{k, 2}, [p.N p.M]);
    Y = X .* H + W;
    [P, rax, vax] = rd_map(p, Y, X, cases{k, 1}, 1, true);
    PdB = 10 * log10(P);
    for mode = {'ca', 'so'}            % the figure shows the last one (SO-CFAR)
        [det, ~] = cfar2d(P, [2 2], [8 6], 1e-6, mode{1});
        % one detection per target: local maxima among detected cells
        Pm = P .* det;
        pk = det & Pm >= imdilate_max(Pm);
        [ir, iv] = find(pk);
        sel = rax(ir) < 120;
        ir = ir(sel);  iv = iv(sel);
        found = zeros(1, numel(tg));
        for t = 1:numel(tg)
            d = abs(rax(ir) - tg(t).R) <= 2 * p.dR & abs(vax(iv) - tg(t).v) <= 2 * p.dv;
            found(t) = any(d);
        end
        ghosts = numel(ir) - sum(found);
        r.rows = [r.rows; k, strcmp(mode{1}, 'so'), found, ghosts, median(PdB(:))];
    end
    ax = subplot(2, 2, k);
    keep = rax <= 100;
    imagesc(vax, rax(keep), PdB(keep, :));  axis xy;  colormap(c.map);  caxis([0 90]);
    hold on;
    for t = 1:numel(tg)
        plot(tg(t).v, tg(t).R, 'o', 'markersize', 16, 'color', c.coral, 'linewidth', 2);
        text(tg(t).v + 4, tg(t).R + 3, tg(t).name, 'color', c.fg, 'fontsize', 10, 'fontweight', 'bold');
    end
    plot(vax(iv), rax(ir), '+', 'color', c.gold, 'markersize', 9, 'linewidth', 1.5);
    xlim([-45 45]);  ylim([0 100]);
    xlabel('radial velocity (m/s)');  ylabel('range (m)');
    ca = r.rows(end-1, :);
    title({tit{k}, sprintf('CA-CFAR %d/5 found, %d false  |  SO-CFAR (+) %d/5, %d false', sum(ca(3:7)), ca(8), sum(found), ghosts)});
    style_axes(ax, c);  set(ax, 'layer', 'top');
    cb = colorbar(ax);  set(cb, 'color', c.muted);  ylabel(cb, 'dB above noise');
    if k == 1, Pdiv = PdB; else, Pcor = PdB; end
end
% range cut through the pedestrian Doppler, and Doppler cut through the car range
ax = subplot(2, 2, 3);  hold on;
[~, iv] = min(abs(vax - tg(3).v));
plot(rax, Pdiv(:, iv), 'color', c.teal, 'linewidth', 1.6);
plot(rax, Pcor(:, iv), 'color', c.coral, 'linewidth', 1.2);
xlim([0 100]);  ylim([-10 95]);
for t = 1:numel(tg), plot([tg(t).R tg(t).R], [-10 95], ':', 'color', c.muted); end
xlabel('range (m)');  ylabel('dB above noise');
title(sprintf('Range cut at v = %.1f m/s (pedestrian bin)', tg(3).v));
style_axes(ax, c);
dark_legend(ax, {'division, 16-QAM', 'correlation, 64-QAM'}, 'northeast');
ax = subplot(2, 2, 4);  hold on;
[~, ir] = min(abs(rax - tg(5).R));
plot(vax, Pdiv(ir, :), 'color', c.teal, 'linewidth', 1.6);
plot(vax, Pcor(ir, :), 'color', c.coral, 'linewidth', 1.2);
xlim([-60 60]);  ylim([-10 95]);
xlabel('radial velocity (m/s)');  ylabel('dB above noise');
title(sprintf('Doppler cut at R = %.0f m (cyclist bin)', tg(5).R));
style_axes(ax, c);
save_figure(h, 'fig1_range_doppler');
r.snr_db = 10 * log10(arrayfun(@(t) abs(t.a)^2, tg));
r.integrated_db = r.snr_db + 10 * log10(p.N * p.M);
write_csv('exp1_targets', 'R_m,v_mps,rcs_m2,snr_per_element_dB,integrated_snr_dB', ...
    [[tg.R].' [tg.v].' [tg.rcs].' r.snr_db.' r.integrated_db.']);
write_csv('exp1_detections', 'case,so_cfar,found_wall,found_car,found_pedestrian,found_truck,found_cyclist,false_alarms,median_floor_dB', r.rows);
end

function D = imdilate_max(P)
% 3x3 maximum filter without the image package
D = P;
for dr = -1:1
    for dc = -1:1
        D = max(D, circshift(P, [dr dc]));
    end
end
end
