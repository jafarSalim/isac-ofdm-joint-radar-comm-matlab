%MAIN  Runs every ISAC experiment and writes figures/ and results/.
%   Works in MATLAB (R2016b or newer, no toolboxes) and GNU Octave (6 or newer).
clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(fullfile(here, 'src'), fullfile(here, 'experiments'));
t0 = tic;
r1 = exp1_scene();
r2 = exp2_crb();
r3 = exp3_modulation();
r4 = exp4_beamforming();
r5 = exp5_detection();
r6 = exp6_numerology();
make_frames();
p = isac_params();
fid = fopen(fullfile(here, 'results', 'summary.txt'), 'w');
fprintf(fid, 'Frame: %d subcarriers x %d symbols, %.1f kHz, B = %.2f MHz, %.2f ms\n', p.N, p.M, p.df/1e3, p.B/1e6, p.M*p.Ts*1e3);
fprintf(fid, 'Resolution %.2f m, %.2f m/s; CP range %.1f m; max velocity %.0f m/s\n', p.dR, p.dv, p.Rcp, p.vmax);
fprintf(fid, 'Scene (found per target, false alarms) [case so wall car ped truck cyclist false floor]:\n%s\n', mat2str(r1.rows, 3));
fprintf(fid, 'Accuracy: RMSE/CRB above -35 dB: range %.2f, velocity %.2f; outliers < 5%% from %g dB; RMSE at -20 dB: %.3g m, %.3g m/s\n', r2.eff, r2.threshold, r2.at20);
fprintf(fid, 'Modulation [floor div, floor corr, PFR div, PFR corr, theory div, theory corr] dB:\n%s\n', mat2str(r3.tab, 3));
fprintf(fid, 'SNR for BER 1e-3 (QPSK, 16, 64): %s dB; rates %s Gb/s\n', mat2str(r3.snr_1e3), mat2str(r3.rate_gbps, 3));
fprintf(fid, 'Beamforming, user rate with the target at -3 dB [60 deg, 5 deg]: Pareto %s, power split %s, weighted sum %s bit/s/Hz; full beam %.2f\n', mat2str(r4.front3, 3), mat2str(r4.split3, 3), mat2str(r4.ws3, 3), r4.ends(1));
fprintf(fid, 'Sensing kept at 90%% rate vs separation %s: %s dB\n', mat2str(r4.seps), mat2str(r4.keep_db, 3));
fprintf(fid, 'Detection: SNR for Pd 0.9 theory/CA/SO %s dB; measured Pfa %s; required SNR (Pd .9, Pfa 1e-6) %.2f dB\n', mat2str(r5.loss, 3), mat2str(r5.pfa_meas, 3), r5.need);
fprintf(fid, 'Detection range [ped cyc car truck] x [20 30 40 dBm]:\n%s\n', mat2str(r5.Rdet, 4));
fprintf(fid, 'CP: peak-to-interference %s dB at %s m\n', mat2str(r6.sinr, 3), mat2str(r6.Rs));
fprintf(fid, 'Run time: %.0f s\n', toc(t0));
fclose(fid);
type(fullfile(here, 'results', 'summary.txt'));
