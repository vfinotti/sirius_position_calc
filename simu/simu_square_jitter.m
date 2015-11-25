clear all
close all
clc

% ---
% System parameters

f_rf    = 499.65794480373e6;    % sirius data frequency
f_adc   = 203/864*f_rf;         % sampling frequency
f_if    = 52/864*f_rf;          % data frequency after sampling

f_rev   = f_rf/864;
f_fofb  = f_rev/6;              % switching frequency;
np      = round(f_adc/f_fofb)*100*3.6;  % number of points in data vector

% ---
% Creating switching signal

t       = (0:np-1)*(1/f_adc);
pm      = sin(2*pi*f_fofb/36*t)*pi/180; % divided by 36 just to have different "f" than "f_fofb"
sw      = square(2*pi*f_fofb*t)*0.25 + 0.75;
sw_j    = square(2*pi*f_fofb*t + pm)*0.25 + 0.75;

% ---
% Carrier signal

sig     = sin(2*pi*f_rf*t);

% ---
% Mixing signals

sig_sw     = sig .* sw;
sig_sw_j   = sig .* sw_j;

% ---
% Plotting 

f = linspace(0,f_adc,np);

figure
plot(t,sw,t,sw_j,'r');
title('Switched Signal - Time Domain')
ylabel('Amplitude')
xlabel('Time (s)')
legend('No Jitter','With Jitter');
axis([0 np/f_adc 0 1.1]);
grid on

semilogy(f,abs(fft(sig_sw)),f,abs(fft(sig_sw_j)),'r');
title('Switched Signal - FFT')
ylabel('FFT')
xlabel('Frequency (Hz)')
legend('No Jitter','With Jitter');
grid on

figure
plot(t,sig_sw,t,sig_sw_j,'r');
title('Jitter Time Domain')
ylabel('Amplitude')
xlabel('Time (s)')
legend('No Jitter','With Jitter');
axis([0 np/f_adc -1.1 1.1]);
grid on

% ---
% Creating Sausaging Window

s           = round(f_adc/(f_fofb*2));    % Number of points on each level of carrier (switching)
s_wind      = zeros(1,s);               % create windows array
zeros_dist  = 5;                        % distance from borders to the transition
size_cos    = 5;                        % size of the cosine, representing the transition
s_half = ceil(s/2);

s_wind(zeros_dist+1:zeros_dist+size_cos)                                = 0.5*cos(linspace(-pi,0,size_cos))+0.5;    % create windows using cos template
s_wind(zeros_dist+size_cos+1:(end-(zeros_dist+size_cos)))               = ones(1,s-2*(size_cos+zeros_dist));        % create pass band
s_wind((end-1)-(zeros_dist+size_cos)+1:(end-1)-zeros_dist)              = 0.5*cos(linspace(0,pi,size_cos))+0.5;     % create windows using cos template

for i = 1:(size(sig_sw_j,2)/s)
    sig_sw_w(s*(i-1)+1:s*(i))   = sig_sw(s*(i-1)+1:s*(i)).*s_wind;        % applying sausaging
    sig_sw_j_w(s*(i-1)+1:s*(i)) = sig_sw_j(s*(i-1)+1:s*(i)).*s_wind;      % applying sausaging
end

%plot(t,mix_j,t,mix_j_w)

% ---
% Cic Filter

s_cos = cos(2*pi*f_fofb*t);
s_sin = sin(2*pi*f_fofb*t);

cic = dfilt.dffir(((1/(s)*ones(1,s)))); % creating fir equivalent to cic

% No jitter - Unsausaged
sig_sw_c     = s_cos.*sig_sw;
sig_sw_s     = -1*s_sin.*sig_sw;
cic_c     = filter(cic,sig_sw_c); % compensated cic
cic_s     = filter(cic,sig_sw_s); % compensated cic
cic_c_d = cic_c(s:s:end);
cic_s_d = cic_s(s:s:end);

% for i = 1:np/s
%     if (mod(i,2) ~= 0 ) % odd number
%         cic_s_d(i)   = cic_s(s*(i-1)+(s_half-1));
%         cic_c_d(i)   = cic_c(s*(i-1)+(s_half-1));
%     else % even number
%         cic_s_d(i)   = cic_s(s*(i-1)+s);
%         cic_c_d(i)   = cic_c(s*(i-1)+s);
%     end
% end

% Jitter - Unsausaged
mix_j_c   = s_cos.*sig_sw_j;
mix_j_s   = -1*s_sin.*sig_sw_j;
cic_j_c   = filter(cic,mix_j_c); % compensated cic
cic_j_s   = filter(cic,mix_j_s); % compensated cic
cic_j_c_d = cic_j_c(s:s:end);
cic_j_s_d = cic_j_s(s:s:end);

% No jitter - Sausaged
mix_w_c     = s_cos.*sig_sw_w;
mix_w_s     = -1*s_sin.*sig_sw_w;
cic_w_c     = filter(cic,mix_w_c); % compensated cic
cic_w_s     = filter(cic,mix_w_s); % compensated cic
cic_w_c_d = cic_w_c(s:s:end);
cic_w_s_d = cic_w_s(s:s:end);

% Jitter - Sausaged
mix_j_w_c   = s_cos.*sig_sw_j_w;
mix_j_w_s   = -1*s_sin.*sig_sw_j_w;
cic_j_w_c   = filter(cic,mix_j_w_c); % compensated cic
cic_j_w_s   = filter(cic,mix_j_w_s); % compensated cic
cic_j_w_c_d = cic_j_w_c(s:s:end);
cic_j_w_s_d = cic_j_w_s(s:s:end);

% ---
% Cordic

% Unsausaged
[cic_j_ph cic_j_amp] = cart2pol(double(cic_j_c_d), double(cic_j_s_d));
[cic_ph cic_amp] = cart2pol(double(cic_c_d), double(cic_s_d));

% Sausaged
[cic_j_w_ph cic_j_w_amp] = cart2pol(double(cic_j_w_c_d), double(cic_j_w_s_d));
[cic_w_ph cic_w_amp] = cart2pol(double(cic_w_c_d), double(cic_w_s_d));

% ---
% Jitter Analysis 

% Unsausaged
[cic_fft, cic_xf] = fourierseries(cic_amp,f_adc/s);
[cic_j_fft, cic_j_xf] = fourierseries(cic_j_amp,f_adc/s);

figure
semilogy(cic_xf, cic_fft, cic_j_xf, cic_j_fft);
title('Unsausaged Signals - FFT')
ylabel('FFT')
xlabel('Frequency (Hz)')
legend('No Jitter', 'with Jitter')
grid on

figure
plot(linspace(0,1/f_fofb*(size(cic_amp,2)),size(cic_amp,2)), cic_amp, linspace(0,1/f_fofb*(size(cic_j_amp,2)),size(cic_j_amp,2)), cic_j_amp);
title('Unsausaged Signals - Amplitude Time Domain')
ylabel('Amplitude')
xlabel('Time (s)')
legend('No jitter', 'with Jitter')
grid on

% Sausaged
[cic_w_fft, cic_w_xf] = fourierseries(cic_w_amp,f_adc/s);
[cic_j_w_fft, cic_j_w_xf] = fourierseries(cic_j_w_amp,f_adc/s);

figure
semilogy(cic_w_xf, cic_w_fft, cic_j_w_xf, cic_j_w_fft);
title('Sausaged Signals - FFT')
ylabel('FFT')
xlabel('Frequency (Hz)')
legend('No Jitter', 'with Jitter')
grid on

figure
plot(linspace(0,1/f_fofb*(size(cic_w_amp,2)),size(cic_w_amp,2)), cic_w_amp, linspace(0,1/f_fofb*(size(cic_j_w_amp,2)),size(cic_j_w_amp,2)), cic_j_w_amp);
title('Sausaged Signals - Amplitude Time Domain')
ylabel('Amplitude')
xlabel('Time (s)')
legend('No Jitter', 'with Jitter')
grid on