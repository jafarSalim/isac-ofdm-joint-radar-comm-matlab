function d = make_frames()
%MAKE_FRAMES  Frames for assets/street.gif: the street seen from the base station.
%   Left: top view (base station at the origin, each road user on its own
%   bearing, moving radially). Right: the range-Doppler map computed from the
%   data frame sent at that instant (division receiver, 16-QAM). Frames are
%   written to <tempdir>/isac_frames and joined into a GIF afterwards.
p = isac_params();  c = palette();
randn('seed', 8);  rand('seed', 8);
d = fullfile(tempdir, 'isac_frames');
if ~exist(d, 'dir'), mkdir(d); end
old = dir(fullfile(d, '*.png'));
for k = 1:numel(old), delete(fullfile(d, old(k).name)); end
brg = [72 8 -38 28 -14];                  % bearings of wall, car, pedestrian, truck, cyclist
nf = 30;  dt = 0.08;
for f = 1:nf
    tg = scene(p, (f - 1) * dt);
    X = qam_mod(16, [p.N p.M]);
    Y = X .* radar_channel(p, tg) + (randn(p.N, p.M) + 1j*randn(p.N, p.M)) / sqrt(2);
    [P, rax, vax] = rd_map(p, Y, X, 'div', 1, true);
    h = figure('visible', 'off', 'position', [0 0 1150 470], 'color', c.bg);
    % top view
    axes('position', [0.03 0.08 0.36 0.84]);  hold on;
    th = linspace(-60, 85, 80);
    for rr = 20:20:100
        plot(rr*cosd(th), rr*sind(th), '-', 'color', c.grid);
    end
    fill([0 100*cosd(th) 0], [0 100*sind(th) 0], c.teal, 'facealpha', 0.06, 'edgecolor', 'none');
    for k = 1:numel(tg)
        x = tg(k).R * cosd(brg(k));  y = tg(k).R * sind(brg(k));
        % radar ping travelling out
        ph = mod(f * 0.35 + k * 0.2, 1);
        plot([0 ph*x], [0 ph*y], '-', 'color', c.teal, 'linewidth', 1);
        plot(x, y, 'o', 'markersize', 6 + 3*log10(tg(k).rcs + 1) * 2, 'color', c.coral, ...
             'markerfacecolor', c.coral);
        text(x + 3, y + 3, sprintf('%s %+.0f m/s', tg(k).name, tg(k).v), 'color', c.fg, 'fontsize', 10);
    end
    % user equipment receiving data
    plot(45*cosd(-62), 45*sind(-62), 's', 'markersize', 10, 'color', c.gold, 'markerfacecolor', c.gold);
    text(45*cosd(-62) + 4, 45*sind(-62), 'phone (16-QAM data)', 'color', c.gold, 'fontsize', 10);
    plot(0, 0, '^', 'markersize', 14, 'color', c.teal, 'markerfacecolor', c.teal);
    axis equal;  xlim([-5 105]);  ylim([-75 105]);  axis off;
    title(sprintf('street, t = %.2f s', (f - 1) * dt), 'color', c.fg, 'fontsize', 13);
    % range-Doppler
    ax = axes('position', [0.47 0.13 0.5 0.76]);
    keep = rax <= 100;
    imagesc(vax, rax(keep), 10*log10(P(keep, :)));  axis xy;  colormap(ax, c.map);  caxis([0 90]);
    xlim([-40 40]);  ylim([0 100]);
    xlabel('radial velocity (m/s)');  ylabel('range (m)');
    title('range-Doppler map from the same OFDM data frame');
    style_axes(ax, c);
    set(h, 'inverthardcopy', 'off', 'PaperPositionMode', 'auto');
    print(h, fullfile(d, sprintf('f%03d.png', f)), '-dpng', '-r90');
    close(h);
end
end
