function a = ula(Nt, theta_deg)
%ULA  Steering vectors of a half-wavelength uniform linear array (Nt x numel(theta)).
a = exp(1j * pi * (0:Nt-1).' * sind(theta_deg(:).'));
end
