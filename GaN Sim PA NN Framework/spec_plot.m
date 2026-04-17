function [Pxx, f] = spec_plot(x, Ts, figNum, title_str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spec_plot (x, Ts, figNum, title_str)
%     Generate power spectrum from voltage waveform data x
%
% Parameters:
%   x          waveform to analyze
%   Ts         timestep for waveform
%   figNum     figure index
%   title_str  figure title
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define sampling frequency
Fs = 1/Ts;  


if length(x) > 4096
    % Define window function
    w = hamming(4096);
else
    % Define window function
    w = hamming(1024);
end


% Compute Welch Power Spectral Density
[Pxx, f] = pwelch(x, w, [], [], Fs, 'centered'); 
f = f / 1e6;  % Convert Hz to MHz

% Normalize power
Pxx = Pxx / max(Pxx);
Pxx = 10 * log10(Pxx);

% Plot Spectrum
figure(figNum);
plot(f, Pxx, 'r');
grid on;
xlabel('Frequency (MHz)');
ylabel('Power Spectral Density (dB)');
title(title_str);

% Define frequency bands for ACPR calculation
main_band = (f >= -35) & (f <= 35);
upper_adj_band = (f > 65) & (f <= 130);
lower_adj_band = (f < -65) & (f >= -130);

% Integrate power over the bands (convert from dB to linear scale)
main_power = sum(10.^(Pxx(main_band) / 10));
upper_adj_power = sum(10.^(Pxx(upper_adj_band) / 10));
lower_adj_power = sum(10.^(Pxx(lower_adj_band) / 10));
adj_power = upper_adj_power + lower_adj_power;

% Calculate ACPR in dB
acpr = 10 * log10(adj_power / main_power);

% Display ACPR result
fprintf('ACPR: %.2f dB\n', acpr);

return;

