function [Pxx, f] = spec_plot (x, Ts, figNum, title_str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%                        
%%%%%%%%                  spec_plot()                   %%%%%%%%
%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spec_plot (x, Ts, figNum)
%     Generate power spectrum from voltage waveform data x
%
%Parameters:
%   x          waveform to analyze
%   Ts         timestep for waveform
%   figNum     figure index
%   title_str  figure title
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paul Draxler
% general utility function
% 3/19/2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(figNum); plot(pwelch(x));

% 
% figure(figNum); plot(periodogram(x));

factor = 1;

if 1
    factor = 1;
    x  = decimate(x, factor, 'fir');
end

Fs = 1/(Ts*factor);

if 0
    XX = fft(x);  XX = fftshift(XX);
    Pxx = 10*log10(XX.* conj(XX) / (length(x)));
else
    w = hamming(4096);
    Pxx = pwelch (x, w, []); 
    Pxx = Pxx/max(abs(Pxx));
    Pxx = 10*log10(fftshift(Pxx));
end

f = Fs*((0:length(Pxx)-1)/length(Pxx) - 0.5)/1e6;

x_max = max(ceil(max(Pxx)/10)*10);

figure(figNum);  
% WiMAX
mask=[-14.75 -50
      -9.75 -32
      -5.45 -25
      -4.75 0
      0 0
      4.75 0
      5.45 -25
      9.75 -32
      14.75 -50];

% Wibro (DL)
mask2=[-14 -60
        -9.27 -60
        -9 -40
        -4.77 -33
        -4.23 0
        0 0 
        4.23 0
        4.77 -33
        9 -40
        9.27 -60
        14 -60];

% Mask for 3GPP LTE 10MHz
mask1=[
      -20 -13-20
      -4 -13-20
      -4 -26-20
      -3.515 -26-20
      -2.715 -14-20
      -2.515 -14-20
      -2.5 0
      0 0
      2.5 0
      2.515 -14-20
      2.715 -14-20
      3.515 -26-20
      4.0 -26-20
      4.0 -13-20
      20 -13-20
      ];

%plot(mask(:,1),mask(:,2),'g',f,Pxx,'r'); grid on
plot(f,Pxx,'r'); grid on
%axis([-15 15 -60 0]);
%plot(f, Pxx);

% if x_max > -10
%     axis([-20 20 -60 x_max]);
% else
%     axis([-20 20 x_max-60 x_max]);
% end

ylabel('Power Spectral Density (dB)');
xlabel('Frequency offset (MHz)');
title(title_str);

Nmin=1;
Nmax=length(x);
Fs_acpr = 1e9;
w_acpr=hann(Nmax-Nmin);
Pxx_acpr=pwelch(x,w_acpr,[]);
Pxx_acpr=fftshift(Pxx_acpr);

f_acpr=linspace(-0.5,0.5,length(Pxx_acpr))*Fs_acpr;

BW_channel = 100e6;
Adj_offset = 120e6;

% Channels index
ind_channel= f_acpr>-1*BW_channel/2 & f_acpr<BW_channel/2;
% ind_upper= f>70e6 & f<200e6;
% ind_lower= f>-200e6 & f<-70e6;

ind_upper= f_acpr>Adj_offset-BW_channel/2 & f_acpr<Adj_offset+BW_channel/2;
ind_lower= f_acpr>-Adj_offset-BW_channel/2 & f_acpr<-Adj_offset+BW_channel/2;

% Riemann integral
int_channel=sum(Pxx_acpr(ind_channel))/sum(ind_channel);
int_upper=sum(Pxx_acpr(ind_upper))/sum(ind_upper);
int_lower=sum(Pxx_acpr(ind_lower))/sum(ind_lower);

acpr_l1=int_lower/int_channel;
acpr_u1=int_upper/int_channel;

acp1_l = 10*log10(acpr_l1);
acp1_u = 10*log10(acpr_u1);
acpr=10*log10(max(acpr_l1,acpr_u1));

fprintf('acpr value: %.2f dB\n', acpr);


return;
