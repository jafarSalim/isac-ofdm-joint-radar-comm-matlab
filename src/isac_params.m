function p = isac_params(mu)
%ISAC_PARAMS  5G NR-like OFDM parameters for joint radar and communication.
%   p = isac_params()   numerology mu = 3 (120 kHz) at 28 GHz, 1024 x 256 frame
%   p = isac_params(mu) any numerology 0..3 with the same 1024 subcarriers.
if nargin < 1, mu = 3; end
p.c      = 299792458;
p.fc     = 28e9;                          % carrier (FR2)
p.lambda = p.c / p.fc;
p.mu     = mu;
p.df     = 15e3 * 2^mu;                   % subcarrier spacing
p.N      = 1024;                          % subcarriers
p.M      = 256;                           % OFDM symbols per frame (2.3 ms)
p.B      = p.N * p.df;                    % occupied bandwidth
p.T      = 1 / p.df;                      % useful symbol time
p.Tcp    = 144 / 2048 * p.T;              % normal cyclic prefix (4.7 % of T)
p.Ts     = p.T + p.Tcp;                   % OFDM symbol with CP
p.Ncp    = round(p.Tcp / p.T * p.N);      % CP samples at fs = N*df
p.dR     = p.c / (2 * p.B);               % range resolution
p.Rmax   = p.c / (2 * p.df);              % unambiguous range (subcarrier spacing)
p.Rcp    = p.c * p.Tcp / 2;               % range covered by the cyclic prefix
p.dv     = p.lambda / (2 * p.M * p.Ts);   % velocity resolution
p.vmax   = p.lambda / (4 * p.Ts);         % unambiguous velocity (+/-)
% link budget of the radar (monostatic, full duplex assumed)
p.Pt_dBm = 30;                            % transmit power
p.G_dBi  = 20;                            % array gain each way
p.NF_dB  = 8;                             % receiver noise figure
p.kT_dBm = -174;                          % dBm/Hz
end
