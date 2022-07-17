function [Y_wave, f] = wave_turb_decomp(u, fs, plot_on)
% decompose the wave velocity from the full velocity signal u. Return
% double-sided Fourier modes of wave velocity signal as complex vector
% containing amplitude and phase info Y_wave = amp_wave*(cos(phs) +
% i*sin(phs)) and frequency f

% interpolate NaNs
u = naninterp(u);

% signal divisions
Ldiv = 256; %size(u,1) - mod(size(u,1),2); %
Ndiv = floor(size(u,1)/Ldiv*2) - 1; %1; %
S2 = zeros(Ldiv,1);
phs2 = zeros(Ldiv,1);

% window and compute fft over each division
hannwin = hann(Ldiv);
for j = 1:Ndiv
    idx_start = (j-1)*Ldiv/2 + 1;
    idx_end = idx_start + Ldiv - 1;
    Y = fft(hannwin.*u( idx_start:idx_end, :));
    S2 = S2 + ((Y).*conj(Y))/(fs*Ldiv); 
    phs2 = phs2 + atan2(imag(Y),real(Y));
end

% average over divisions and convert to single-sided spectrum
S2 = S2/Ndiv; %two-sided spectrum
% S2 = [S2(Ldiv/2+1:Ldiv); S2(1:Ldiv/2+1)];
S1 = S2(1:Ldiv/2+1); S1(1:end) = S1(1:end)*2; % single sided spectrum
phs2 = phs2/Ndiv;
phs2 = [phs2(Ldiv/2+1:Ldiv); phs2(1:Ldiv/2+1)];
% phs1 = phs2(1:Ldiv/2+1);

% frequency range
f1 = fs*(0:(Ldiv/2))/Ldiv; % frequncy range for single-sided spectrum
f = fs*(-(Ldiv/2):(Ldiv/2))'/Ldiv; % frequency vector for two-sided spectrum
% plot to see wave peak
if plot_on
    figure; loglog(f,S2,'b-'); xlabel('f [Hz]'); ylabel('E'); 
    [~,A_idx] = min(abs(f-1)); A = S2(A_idx);
    hold on; loglog(f(f>.1),A*f(f>.1).^(-5/3),'-k')
end

% l and r bounds of wave peak [Hz]
wavepeak_l = 1.2;
wavepeak_r = 3;
wave_idx = f>wavepeak_l & f<wavepeak_r;
wave_idx_neg = f<-wavepeak_l & f>-wavepeak_r;  % negative side of spectrum

% l and r bound of turbulence range (inertial subrange) [Hz]
turbrange_l = 0.5;%4.5;
turbrange_r = 20;
turb_idx = f>turbrange_l & f<turbrange_r;

% interpolate spectrum over wave peak
log_f_turb = log(f(turb_idx & ~wave_idx));
log_S2_turb = log(S2(turb_idx & ~wave_idx));
P = polyfit(log_f_turb,log_S2_turb,1);
S2_turb_hat = exp(polyval(P,log(abs(f))));

if plot_on
    hold on; plot(f,S2_turb_hat,'r-');
    figure; plot(f,S2,'b-',f,S2_turb_hat,'r-');
end

% wave spectrum
S2_wave = zeros(size(S2));
S2_wave(wave_idx | wave_idx_neg) = S2(wave_idx | wave_idx_neg) - S2_turb_hat(wave_idx | wave_idx_neg);

% recover Fourier coeff amplitudes
% S2_wave = [S1_wave(Ldiv:-1:2); S1_wave];
amp2_wave = sqrt(S2_wave); %*Ldiv;

Y_wave = amp2_wave.*(cos(phs2) + 1j*sin(phs2));
% if I change Ndiv to 1, I can recover the full wave velocity time series
% as u_wave = ifft(Y_wave), then take the hilbert transform of it to get
% wave phase

end
function X = naninterp(X)
% Interpolate over NaNs
% See INTERP1 for more info
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)));
end