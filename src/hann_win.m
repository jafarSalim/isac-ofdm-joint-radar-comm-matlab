function w = hann_win(n)
%HANN_WIN  Symmetric Hann window (column), no toolbox needed.
w = 0.5 - 0.5 * cos(2 * pi * (0:n-1).' / (n - 1));
end
