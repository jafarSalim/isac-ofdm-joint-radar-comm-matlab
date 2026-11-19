function c = palette()
%PALETTE  "Reef" colours: teal sensing, coral communication, on charcoal.
c.bg     = [0.06 0.08 0.10];
c.panel  = [0.09 0.12 0.15];
c.fg     = [0.90 0.95 0.96];
c.grid   = [0.25 0.32 0.36];
c.muted  = [0.55 0.63 0.67];
c.teal   = [0.00 0.83 0.75];
c.coral  = [1.00 0.42 0.38];
c.gold   = [1.00 0.80 0.20];
c.violet = [0.66 0.48 1.00];
c.sky    = [0.30 0.67 1.00];
c.lime   = [0.60 0.90 0.25];
c.mods   = {c.teal, c.gold, c.coral};          % QPSK, 16QAM, 64QAM
c.lw     = 2.4;
k = [0 0.03 0.04 0.07; 0.25 0.00 0.25 0.35; 0.5 0.00 0.60 0.62; ...
     0.72 0.55 0.85 0.55; 0.88 1.00 0.75 0.40; 1 1.00 0.98 0.90];
x = linspace(0, 1, 256)';
c.map = [interp1(k(:,1),k(:,2),x) interp1(k(:,1),k(:,3),x) interp1(k(:,1),k(:,4),x)];
end
