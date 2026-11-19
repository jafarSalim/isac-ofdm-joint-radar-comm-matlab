function bits = qam_demod(M, y)
%QAM_DEMOD  Hard-decision Gray demapper matching QAM_MOD.
y = y(:);  k = log2(M);
if M == 2
    bits = real(y) > 0;
    return;
end
L = sqrt(M);  h = k / 2;
s = y * sqrt(2 * (M - 1) / 3);
bI = lvl2gray(real(s), L, h);  bQ = lvl2gray(imag(s), L, h);
bits = reshape([bI bQ].', [], 1);
end

function b = lvl2gray(v, L, h)
idx = round((v + (L - 1)) / 2);
idx = min(max(idx, 0), L - 1);
b = false(numel(v), h);
for i = h:-1:1
    b(:, i) = mod(idx, 2);  idx = floor(idx / 2);
end
% binary -> Gray
for i = h:-1:2
    b(:, i) = xor(b(:, i), b(:, i-1));
end
end
