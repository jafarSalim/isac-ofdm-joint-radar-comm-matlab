function [x, bits] = qam_mod(M, sz, bits)
%QAM_MOD  Gray-coded square M-QAM with unit average energy.
%   [x, bits] = qam_mod(M, [rows cols])          random data
%   x = qam_mod(M, [rows cols], bits)             given bits (prod(sz)*log2(M) x 1)
k = log2(M);  n = prod(sz);
if nargin < 3, bits = rand(n * k, 1) > 0.5; end
b = reshape(bits, k, n).';
if M == 2
    s = 2 * b - 1;
    x = reshape(s, sz);
    return;
end
h = k / 2;  L = sqrt(M);
I = gray2lvl(b(:, 1:h), L);  Q = gray2lvl(b(:, h+1:end), L);
s = (I + 1j * Q) / sqrt(2 * (M - 1) / 3);
x = reshape(s, sz);
end

function lvl = gray2lvl(b, L)
% Gray bits -> PAM level in {-(L-1), ..., L-1}
g = b(:, 1);  idx = g;
for i = 2:size(b, 2)
    g = xor(g, b(:, i));
    idx = 2 * idx + g;
end
lvl = 2 * double(idx) - (L - 1);
end
