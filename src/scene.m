function tg = scene(p, t)
%SCENE  Street scene seen by a roadside 28 GHz base station at time t (s).
%   Returns targets with range R (m), radial velocity v (m/s, + = receding),
%   radar cross-section rcs (m^2), echo amplitude a and a name.
if nargin < 2, t = 0; end
names = {'glass wall', 'car', 'pedestrian', 'truck', 'cyclist'};
R0    = [14     38     48     72    80];
v     = [0     -18    1.4     12     6];
rcs   = [300    10    0.5    100     1];
ph    = [0.3    1.9    4.1    2.2   5.0];          % fixed phases (reproducible)
tg = struct('name', {}, 'R', {}, 'v', {}, 'rcs', {}, 'a', {});
for k = 1:numel(R0)
    R = R0(k) + v(k) * t;
    snr = echo_snr(p, R, rcs(k));
    tg(k).name = names{k};  tg(k).R = R;  tg(k).v = v(k);  tg(k).rcs = rcs(k);
    tg(k).a = sqrt(snr) * exp(1j * ph(k));          % noise has unit variance
end
end
