clear all
%close all
%clc

%% Parameters Settings

% ------------------------------------------------------------
% Signal

f_rf    = 499.65794480373e6;
f_adc   = 203/864*f_rf;                 % sampling frequency
f_if    = 52/864*f_rf;                  % data frequency after sampling
f_rev   = f_rf/864;
s_f_fofb  = f_rev/6;
s_f_sw    = s_f_fofb;                       % switching frequency;
s_np_sw   = round(f_adc/s_f_sw);
s_n_sw    = 1000;                         % number of switching periods
s_np      = s_np_sw*s_n_sw;                   % number of points in data vector

% T       = [1.0 1.35 1.7 1.4];         % gains  T1, T2, T3 and T4 (Pairs are T1/T3, T2/T4)
T       = [1.1028 0.892296 1.0631 0.955918];         % gains  T1, T2, T3 and T4 (Pairs are T1/T3, T2/T4)
% T       = [1 1 1 1];                    % gains  T1, T2, T3 and T4 (Pairs are T1/T3, T2/T4)
phases  = [-1.21 -0.70 -1.17 -1.00]* 0.5;    % phase delays in rad
% phases  = [-pi -0 -pi -0];    % phase delays in rad
s_size_exp = 200;                           % size of the exponential, representing the switching transition

noise_flag = 0;                         % 1 for "noise", 0 for "without noise"

beam_xy = [-1 1]*1e3;                  % position of the beam in  nm, [x y]
button_r = 3e6;                         % button radius in nm
chamber_r = 12e6;                       % chamber radius in nm
Kx = 8.5785e6;                          % converts to nm

% ------------------------------------------------------------
% Sirius

s_filter_flag = 'fir';                  % 'fir' or 'cic'

% ------------------------------------------------------------
% Libera

l_f_fofb  = f_rev/43;
l_f_sw    = l_f_fofb;                   % switching frequency;
l_np_sw   = round(f_adc/l_f_sw);
l_n_sw    = 100;                         % number of switching periods
l_np      = l_np_sw*l_n_sw*4;               % number of points in data vector
l_size_exp = 0;%round(l_np_sw/1000);                           % size of the exponential, representing the switching transition

l_ph_comp = 1;
l_amp_comp = 1;
l_spk_rm = 0;
l_filter_flag = 'fir';                  % 'fir' or 'firn'


%% Creating signals - Sirius

% ---
% Creating time signal and frequency array using the sampling frequency

s_t       = (0:s_np-1)*(1/f_adc);
s_fx      = linspace(0,f_adc,s_np);

% Line gains specified on the beggining

% ---
% Data signal

abcd = pos2abcd(beam_xy,button_r,chamber_r); % obtaining the antennas intensity for the considered position

% phase delays specified in on the beginning

for i = 0:s_np/s_np_sw*2-1
    if rem(i,2) == 0 % even, delay of its channel
        s_sig(i*s_np_sw/2+1:((i+1)*s_np_sw)/2,1) = abcd(1) * sin(2*pi*f_rf*s_t(i*s_np_sw/2+1:((i+1)*s_np_sw)/2)+phases(1)); % Signal for line 1
        s_sig(i*s_np_sw/2+1:((i+1)*s_np_sw)/2,2) = abcd(2) * sin(2*pi*f_rf*s_t(i*s_np_sw/2+1:((i+1)*s_np_sw)/2)+phases(2)); % Signal for line 2
        s_sig(i*s_np_sw/2+1:((i+1)*s_np_sw)/2,3) = abcd(3) * sin(2*pi*f_rf*s_t(i*s_np_sw/2+1:((i+1)*s_np_sw)/2)+phases(3)); % Signal for line 3
        s_sig(i*s_np_sw/2+1:((i+1)*s_np_sw)/2,4) = abcd(4) * sin(2*pi*f_rf*s_t(i*s_np_sw/2+1:((i+1)*s_np_sw)/2)+phases(4)); % Signal for line 4
    else % odd, delay of opposing channel
        s_sig(i*s_np_sw/2+1:((i+1)*s_np_sw)/2,1) = abcd(1) * sin(2*pi*f_rf*s_t(i*s_np_sw/2+1:((i+1)*s_np_sw)/2)+phases(3)); % Signal for line 1
        s_sig(i*s_np_sw/2+1:((i+1)*s_np_sw)/2,2) = abcd(2) * sin(2*pi*f_rf*s_t(i*s_np_sw/2+1:((i+1)*s_np_sw)/2)+phases(4)); % Signal for line 2
        s_sig(i*s_np_sw/2+1:((i+1)*s_np_sw)/2,3) = abcd(3) * sin(2*pi*f_rf*s_t(i*s_np_sw/2+1:((i+1)*s_np_sw)/2)+phases(1)); % Signal for line 3
        s_sig(i*s_np_sw/2+1:((i+1)*s_np_sw)/2,4) = abcd(4) * sin(2*pi*f_rf*s_t(i*s_np_sw/2+1:((i+1)*s_np_sw)/2)+phases(2)); % Signal for line 4
    end
end

% ---
% Adding noise
if noise_flag
    noise   = randn(s_np,4)*1e-3;
    s_sig     = s_sig + noise;
end

% ---
% Creating Unswitched signal

s_sw0(:,1)   = square(2*pi*s_f_fofb*s_t(1:s_np_sw))*(T(1)-T(3))/2 + T(3) +(T(1)-T(3))/2;           % Signal for line 1
s_sw0(1:s_size_exp,1)   = (T(3)-T(1))*exp(-linspace(0,6,s_size_exp)) + T(1);                       % Transition 1
s_sw0(s_np_sw/2+1:s_np_sw/2+s_size_exp,1)   = (T(1)-T(3))*exp(-linspace(0,6,s_size_exp)) + T(3);   % Transition 2

s_sw0(:,2)   = square(2*pi*s_f_fofb*s_t(1:s_np_sw))*(T(2)-T(4))/2 + T(4) +(T(2)-T(4))/2;           % Signal for line 2
s_sw0(1:s_size_exp,2)   = (T(4)-T(2))*exp(-linspace(0,6,s_size_exp)) + T(2);                       % Transition 1
s_sw0(s_np_sw/2+1:s_np_sw/2+s_size_exp,2)   = (T(2)-T(4))*exp(-linspace(0,6,s_size_exp)) + T(4);   % Transition 2

s_sw0(:,3)   = square(2*pi*s_f_fofb*s_t(1:s_np_sw)+pi)*(T(1)-T(3))/2 + T(3) +(T(1)-T(3))/2;        % Signal for line 3
s_sw0(1:s_size_exp,3)   = (T(1)-T(3))*exp(-linspace(0,6,s_size_exp)) + T(3);                       % Transition 1
s_sw0(s_np_sw/2+1:s_np_sw/2+s_size_exp,3)   = (T(3)-T(1))*exp(-linspace(0,6,s_size_exp)) + T(1);   % Transition 2

s_sw0(:,4)   = square(2*pi*s_f_fofb*s_t(1:s_np_sw)+pi)*(T(2)-T(4))/2 + T(4) +(T(2)-T(4))/2;        % Signal for line 4
s_sw0(1:s_size_exp,4)   = (T(2)-T(4))*exp(-linspace(0,6,s_size_exp)) + T(4);                       % Transition 1
s_sw0(s_np_sw/2+1:s_np_sw/2+s_size_exp,4)   = (T(4)-T(2))*exp(-linspace(0,6,s_size_exp)) + T(2);   % Transition 2

s_sw = repmat(s_sw0,s_n_sw,1);

s_sig_sw = zeros(s_np_sw*s_n_sw,4); % allocating variable
for i = 1:size(s_sig,2)
    s_sig_sw(:,i)  = s_sig(:,i) .* s_sw(:,i);
end

% ---
% Plotting switched signal
figure
plot(1:s_np_sw,s_sw0(:,1),1:s_np_sw,square(2*pi*s_f_fofb*s_t(1:s_np_sw))*(T(1)-T(3))/2 + T(3) +(T(1)-T(3))/2,'r');
title('Switched Signal (Sirius) - Channel 1 (A and C) Amplitudes');
ylabel('Amplitude')
xlabel('Time')
grid on
legend('Attenuated','Ideal','location','best')
axis([1 size(s_sw0,1) min(s_sw0(:,1))*0.95 max(s_sw0(:,1))*1.05]);

% ---
% Plotting switched signal
figure
plot(s_t,s_sig_sw(:,1),s_t,s_sw(:,1),'r');
title('Unswitched Signal (Sirius) - Channel 1 (A and C)');
ylabel('Amplitude')
xlabel('Time [s]')
legend('Unswitched Signal','Switching','location','best')
grid on

%% Creating signals - Libera

% ---
% Creating time signal and frequency array using the sampling frequency

l_t       = (0:l_np-1)*(1/f_adc);
l_fx      = linspace(0,f_adc,l_np);

% Line gains specified on the beggining

% ---
% Data signal

abcd = pos2abcd(beam_xy,button_r,chamber_r); % obtaining the antennas intensity for the considered position

% phase delays specified in on the beginning
l_sig = zeros(l_np,4);
if l_ph_comp
%     l_sig(:,1) = abcd(1) * sin(2*pi*f_rf*l_t + phases(1)); % Signal for line 1
%     l_sig(:,2) = abcd(2) * sin(2*pi*f_rf*l_t + phases(1)); % Signal for line 2
%     l_sig(:,3) = abcd(3) * sin(2*pi*f_rf*l_t + phases(1)); % Signal for line 3
%     l_sig(:,4) = abcd(4) * sin(2*pi*f_rf*l_t + phases(1)); % Signal for line 4
    l_sig(:,1) = abcd(1) * sin(2*pi*f_rf*l_t + phases(1)); % Signal for line 1
    l_sig(:,2) = abcd(2)/abcd(1) * l_sig(:,1); % Signal for line 2
    l_sig(:,3) = abcd(3)/abcd(1) * l_sig(:,1); % Signal for line 3
    l_sig(:,4) = abcd(4)/abcd(1) * l_sig(:,1); % Signal for line 4

else
%     l_sig0(:,1) = sin(2*pi*f_rf*l_t + phases(1)); % Signal for line 1
%     l_sig0(:,2) = sin(2*pi*f_rf*l_t + phases(1)); % Signal for line 2
%     l_sig0(:,3) = sin(2*pi*f_rf*l_t + phases(1)); % Signal for line 3
%     l_sig0(:,4) = sin(2*pi*f_rf*l_t + phases(1)); % Signal for line 4
    
    l_sig0(1:l_np_sw,1) = abcd(1)*sin(2*pi*f_rf*l_t(1:l_np_sw) + phases(3)); % Switching Stage 1
    l_sig0(1:l_np_sw,2) = abcd(2)*sin(2*pi*f_rf*l_t(1:l_np_sw) + phases(4));
    l_sig0(1:l_np_sw,3) = abcd(3)*sin(2*pi*f_rf*l_t(1:l_np_sw) + phases(1));
    l_sig0(1:l_np_sw,4) = abcd(4)*sin(2*pi*f_rf*l_t(1:l_np_sw) + phases(2));

    l_sig0(l_np_sw+1:2*l_np_sw,1) = abcd(1)*sin(2*pi*f_rf*l_t(l_np_sw+1:2*l_np_sw) + phases(2)); % Switching Stage 2
    l_sig0(l_np_sw+1:2*l_np_sw,2) = abcd(2)*sin(2*pi*f_rf*l_t(l_np_sw+1:2*l_np_sw) + phases(1));
    l_sig0(l_np_sw+1:2*l_np_sw,3) = abcd(3)*sin(2*pi*f_rf*l_t(l_np_sw+1:2*l_np_sw) + phases(4));
    l_sig0(l_np_sw+1:2*l_np_sw,4) = abcd(4)*sin(2*pi*f_rf*l_t(l_np_sw+1:2*l_np_sw) + phases(3));

    l_sig0(2*l_np_sw+1:3*l_np_sw,1) = abcd(1)*sin(2*pi*f_rf*l_t(2*l_np_sw+1:3*l_np_sw) + phases(4)); % Switching Stage 3
    l_sig0(2*l_np_sw+1:3*l_np_sw,2) = abcd(2)*sin(2*pi*f_rf*l_t(2*l_np_sw+1:3*l_np_sw) + phases(3));
    l_sig0(2*l_np_sw+1:3*l_np_sw,3) = abcd(3)*sin(2*pi*f_rf*l_t(2*l_np_sw+1:3*l_np_sw) + phases(2));
    l_sig0(2*l_np_sw+1:3*l_np_sw,4) = abcd(4)*sin(2*pi*f_rf*l_t(2*l_np_sw+1:3*l_np_sw) + phases(1));

    l_sig0(3*l_np_sw+1:4*l_np_sw,1) = abcd(1)*sin(2*pi*f_rf*l_t(3*l_np_sw+1:4*l_np_sw) + phases(1)); % Switching Stage 4
    l_sig0(3*l_np_sw+1:4*l_np_sw,2) = abcd(2)*sin(2*pi*f_rf*l_t(3*l_np_sw+1:4*l_np_sw) + phases(2));
    l_sig0(3*l_np_sw+1:4*l_np_sw,3) = abcd(3)*sin(2*pi*f_rf*l_t(3*l_np_sw+1:4*l_np_sw) + phases(3));
    l_sig0(3*l_np_sw+1:4*l_np_sw,4) = abcd(4)*sin(2*pi*f_rf*l_t(3*l_np_sw+1:4*l_np_sw) + phases(4));

    l_sig = repmat(l_sig0,l_n_sw,1);
end
% ---
% Adding noise
if noise_flag
    noise   = randn(l_np,4)*1e-3;
    l_sig     = l_sig + noise;
end

% ---
% Creating Unswitched signal

l_sw0(1:2*l_np_sw,1)   = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(3)-T(2))/2 + T(2) +(T(3)-T(2))/2;                % Signal for line 1
l_sw0(2*l_np_sw+1:4*l_np_sw,1)   = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(4)-T(1))/2 + T(1) +(T(4)-T(1))/2;      % Signal for line 1
l_sw0(1:l_size_exp,1)   = (T(1)-T(3))*exp(-linspace(0,6,l_size_exp)) + T(3);                        % Transition 1
l_sw0(l_np_sw+1:l_np_sw+l_size_exp,1)   = (T(3)-T(2))*exp(-linspace(0,6,l_size_exp)) + T(2);        % Transition 2
l_sw0(2*l_np_sw+1:2*l_np_sw+l_size_exp,1)   = (T(2)-T(4))*exp(-linspace(0,6,l_size_exp)) + T(4);    % Transition 3
l_sw0(3*l_np_sw+1:3*l_np_sw+l_size_exp,1)   = (T(4)-T(1))*exp(-linspace(0,6,l_size_exp)) + T(1);    % Transition 4

l_sw0(1:2*l_np_sw,2)   = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(4)-T(1))/2 + T(1) +(T(4)-T(1))/2;                % Signal for line 2
l_sw0(2*l_np_sw+1:4*l_np_sw,2)   = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(3)-T(2))/2 + T(2) +(T(3)-T(2))/2;      % Signal for line 2
l_sw0(1:l_size_exp,2)   = (T(2)-T(4))*exp(-linspace(0,6,l_size_exp)) + T(4);                        % Transition 1
l_sw0(l_np_sw+1:l_np_sw+l_size_exp,2)   = (T(4)-T(1))*exp(-linspace(0,6,l_size_exp)) + T(1);        % Transition 2
l_sw0(2*l_np_sw+1:2*l_np_sw+l_size_exp,2)   = (T(1)-T(3))*exp(-linspace(0,6,l_size_exp)) + T(3);    % Transition 3
l_sw0(3*l_np_sw+1:3*l_np_sw+l_size_exp,2)   = (T(3)-T(2))*exp(-linspace(0,6,l_size_exp)) + T(2);    % Transition 4

l_sw0(1:2*l_np_sw,3)   = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(1)-T(4))/2 + T(4) +(T(1)-T(4))/2;                % Signal for line 3
l_sw0(2*l_np_sw+1:4*l_np_sw,3)   = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(2)-T(3))/2 + T(3) +(T(2)-T(3))/2;      % Signal for line 3
l_sw0(1:l_size_exp,3)   = (T(3)-T(1))*exp(-linspace(0,6,l_size_exp)) + T(1);                        % Transition 1
l_sw0(l_np_sw+1:l_np_sw+l_size_exp,3)   = (T(1)-T(4))*exp(-linspace(0,6,l_size_exp)) + T(4);        % Transition 2
l_sw0(2*l_np_sw+1:2*l_np_sw+l_size_exp,3)   = (T(4)-T(2))*exp(-linspace(0,6,l_size_exp)) + T(2);    % Transition 3
l_sw0(3*l_np_sw+1:3*l_np_sw+l_size_exp,3)   = (T(2)-T(3))*exp(-linspace(0,6,l_size_exp)) + T(3);    % Transition 4

l_sw0(1:2*l_np_sw,4)   = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(2)-T(3))/2 + T(3) +(T(2)-T(3))/2;                % Signal for line 4
l_sw0(2*l_np_sw+1:4*l_np_sw,4)   = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(1)-T(4))/2 + T(4) +(T(1)-T(4))/2;      % Signal for line 4
l_sw0(1:l_size_exp,4)   = (T(2)-T(2))*exp(-linspace(0,6,l_size_exp)) + T(2);                        % Transition 1
l_sw0(l_np_sw+1:l_np_sw+l_size_exp,4)   = (T(2)-T(3))*exp(-linspace(0,6,l_size_exp)) + T(3);        % Transition 2
l_sw0(2*l_np_sw+1:2*l_np_sw+l_size_exp,4)   = (T(3)-T(1))*exp(-linspace(0,6,l_size_exp)) + T(1);    % Transition 3
l_sw0(3*l_np_sw+1:3*l_np_sw+l_size_exp,4)   = (T(1)-T(4))*exp(-linspace(0,6,l_size_exp)) + T(4);    % Transition 4

l_sw = repmat(l_sw0,l_n_sw,1);

l_sig_sw = zeros(l_np_sw*l_n_sw*4,4); % allocating variable
for i = 1:size(l_sig,2)
    l_sig_sw(:,i)  = l_sig(:,i) .* l_sw(:,i);
end

% ---
% Plotting switched signal
s_temp1 = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(3)-T(2))/2 + T(2) +(T(3)-T(2))/2;
s_temp2 = square(2*pi*l_f_sw/2*l_t(1:2*l_np_sw))*(T(4)-T(1))/2 + T(1) +(T(4)-T(1))/2;

figure
plot(1:4*l_np_sw,l_sw0(:,1),1:4*l_np_sw,[s_temp1 s_temp2],'r');
title('Switched Signal (Libera) - Channel 1 (A and C) Amplitudes');
ylabel('Amplitude')
xlabel('Time')
grid on
legend('Attenuated','Ideal','location','best')
axis([1 size(l_sw0,1) min(l_sw0(:,1))*0.95 max(l_sw0(:,1))*1.05]);

figure
semilogy(s_fx,abs(fft(s_sig_sw(:,1))),l_fx,abs(fft(l_sig_sw(:,1))),'r')
title('Unswitched Signal - Channel 1 (A and C)');
ylabel('FFT')
xlabel('Frequency [Hz]')
legend('Sirius','Libera','location','best')
grid on

%% Processing - Sirius

% s_mix_c = zeros(np_sw*n_sw,4);    % allocating variable
% s_mix_s = zeros(np_sw*n_sw,4);    % allocating variable
%
% s_cic_c = zeros(np_sw*n_sw,4);    % allocating variable
% s_cic_s = zeros(np_sw*n_sw,4);    % allocating variable
%
% s_cic_c_d = zeros(2*n_sw,4);        % allocating variable
% s_cic_s_d = zeros(2*n_sw,4);        % allocating variable
%
% s_cic_ph0 = zeros(2*n_sw,4);        % allocating variable
% s_cic_amp0 = zeros(2*n_sw,4);       % allocating variable
%
% s_cic_ph = zeros(2*n_sw-1,4);       % allocating variable
% s_cic_amp = zeros(2*n_sw-1,4);      % allocating variable
%
% s_cic_amp_am = zeros(2*n_sw-1,4);   % allocating variable

for i = 1:size(s_sig,2)

    % ---
    % Mixing signal
    s_cos       = cos(2*pi*f_if*s_t);
    s_sin       = sin(2*pi*f_if*s_t);
    s_mix_c(:,i)  = s_sig_sw(:,i).*s_cos';
    s_mix_s(:,i)  = -1*s_sig_sw(:,i).*s_sin';
end

% ---
% Filter
s_wind_size   = round(f_adc/(s_f_sw));      % size of the windows

s_cic             = dfilt.dffir(((1/(s_wind_size)*ones(1,s_wind_size)))); % creating fir equivalent to cic

for i = 1:size(s_sig,2)
    if strcmp(s_filter_flag ,'fir')
        % Fir
        s_cic_c(:,i)      = filter(s_cic,s_mix_c(:,i)); % cic
        s_cic_s(:,i)      = filter(s_cic,s_mix_s(:,i)); % cic
        s_cic_c_d(:,i)    = s_cic_c(s_wind_size:s_wind_size:end,i); % decimating
        s_cic_s_d(:,i)    = s_cic_s(s_wind_size:s_wind_size:end,i); % decimating

    elseif strcmp(s_filter_flag ,'cic')
        % CIC
        s_m = 1;            % Differential delays in the filter.
        s_n = 1;            % Filter sections
        s_r = s_wind_size;  % Decimation factor
        s_hm          = mfilt.cicdecim(s_r,s_m,s_n);             % Expects 16-bit input by default.

        s_fp_max = 1.125;   % max number in float point conversion
        s_cic_s_d_32  = filter(s_hm, int32(s_mix_s*(2^16)/s_fp_max));
        s_cic_c_d_32  = filter(s_hm, int32(s_mix_c*(2^16)/s_fp_max));
        s_cic_s_d     = double(s_cic_s_d_32)*s_fp_max/(2^16);  % converts back to float
        s_cic_c_d     = double(s_cic_c_d_32)*s_fp_max/(2^16);  % converts back to float
    end
end

for i = 1:size(s_sig,2)
    % ---
    % Cordic
    [s_cic_ph0(:,i), s_cic_amp0(:,i)] = cart2pol(s_cic_c_d(:,i), s_cic_s_d(:,i));

    s_cic_amp(:,i) = s_cic_amp0(2:end,i); % correcting cic beginning
    s_cic_ph(:,i) = s_cic_ph0(2:end,i); % correcting cic beginning

    % Filtered Parameters
    s_f_cic   = round(f_adc/(s_wind_size));                    % updating decimated frequency
    s_t_cic   = linspace(0,max(s_t),length(s_cic_amp));

end

% ---
% Calculating X and Y - D/S and P D/S

[s_cic_xy_pds]        = calcpos_pds(s_cic_amp(1:end,:),Kx);       % PDS
[s_cic_xy_ds]         = calcpos_ds(s_cic_amp(1:end,:),Kx);        % DS

% picking "x" values
s_cic_x_pds = s_cic_xy_pds(:,1) - mean(s_cic_xy_pds(:,1));
s_cic_x_ds = s_cic_xy_ds(:,1) - mean(s_cic_xy_ds(:,1));



% Plotting

% ---
% Comparing integrated noise

% PSD calculation
s_wind    = rectwin(s_n_sw/2);
s_nover   = s_n_sw/4;
s_w       = s_n_sw/2;

s_cic_x_pds_psd       = pwelch(s_cic_x_pds,s_wind,s_nover,s_w,s_f_cic);
s_cic_x_ds_psd       = pwelch(s_cic_x_ds,s_wind,s_nover,s_w,s_f_cic);

% Integrating signals

s_f_step = s_f_cic/(2*length(s_cic_x_pds_psd)); % frequency/Hz step

s_cic_x_ds_psd_int       = cumtrapz(s_cic_x_ds_psd)*s_f_step;
s_cic_x_pds_psd_int       = cumtrapz(s_cic_x_pds_psd)*s_f_step;

%% Processing - Libera

if l_amp_comp
    if l_ph_comp
        for i = 0:3
            if i == 0 % switching position 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(3)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(4)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(1)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(2)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(3)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(4)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(1)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(2)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 4
                
            elseif i== 1 % switching position 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(2)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(1)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(4)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(3)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(2)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(1)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(4)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(3)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 4
            elseif i== 2 % switching position 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(4)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(3)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(2)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(1)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(4)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(3)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(2)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(1)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 4
            else % switching position 4
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(1)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(2)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(3)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(4)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(1)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(2)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(3)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(4)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1))); % Signal for line 4
            end
        end
    else
        for i = 0:3
            if i == 0 % switching position 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(3)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(4)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(1)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(2)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(3)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(4)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(1)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(2)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 4
                
            elseif i== 1 % switching position 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(2)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(1)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(4)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(3)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(2)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(1)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(4)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(3)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 4
            elseif i== 2 % switching position 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(4)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(3)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(2)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(1)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(4)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(3)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(2)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(1)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 4
            else % switching position 4
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(1)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(2)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(3)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(4)*cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = T(1)/T(1)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = T(1)/T(2)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = T(1)/T(3)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = T(1)/T(4)*sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw))); % Signal for line 4
            end
        end
    end
else
    if l_ph_comp
        for i = 0:3
            if i == 0 % switching position 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 4
                
            elseif i== 1 % switching position 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 4
            elseif i== 2 % switching position 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 4
            else % switching position 4
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) + (phases(1)) ); % Signal for line 4
            end
            
        end
    else
        
        for i = 0:3
            if i == 0 % switching position 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 4
                
            elseif i== 1 % switching position 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 4
            elseif i== 2 % switching position 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 4
            else % switching position 4
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),1) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 1
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),2) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 2
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),3) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 3
                l_cos0(i*l_np_sw+1:((i+1)*l_np_sw),4) = cos(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 4
                
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),1) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 1
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),2) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 2
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),3) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 3
                l_sin0(i*l_np_sw+1:((i+1)*l_np_sw),4) = sin(2*pi*f_if*l_t(i*l_np_sw+1:((i+1)*l_np_sw)) ); % Signal for line 4
            end
        end
        
    end
end

l_sin = repmat(l_sin0,l_n_sw,1); % repeating l_sin0 and l_cos0
l_cos = repmat(l_cos0,l_n_sw,1); % repeating l_sin0 and l_cos0

l_mix_c = zeros(size(l_sig_sw,1),4);
l_mix_s = l_mix_c;

for i = 1:size(l_sig,2)

    % ---
    % Mixing signal & amplitude and phase compensation
    l_mix_c(:,i)  = l_sig_sw(:,i).*l_cos(:,i);
    l_mix_s(:,i)  = -1*l_sig_sw(:,i).*l_sin(:,i);
end

% ---
% Filter
l_wind_size   = round(l_np_sw*4/3); %round(f_adc/(l_f_sw)/l_f_sw*10000); %round(f_adc/(l_f_sw));      % size of the windows - number of points calculated to achieve 10 KHz

% l_fir             = dfilt.dffir(((1/(l_wind_size)*ones(1,l_wind_size)))); % creating fir equivalent to cic

%  l_fir             = ((1/(l_wind_size)*ones(1,l_wind_size))); % creating fir equivalent to cic
load l_fir;
load l_dec; % Nyquist decimator filter of "l_wind_size" decimation factor

l_fir_c = zeros(size(l_mix_c,1),4);
l_fir_s = l_fir_c;
l_fir_c_d = zeros(floor(size(l_mix_c,1)/l_wind_size),4);
l_fir_s_d = l_fir_c_d;
l_fir_s_d_n = l_fir_c_d;
l_fir_c_d_n = l_fir_c_d;
l_fir_ph0 = l_fir_c_d;
l_fir_amp0 = l_fir_c_d;
l_fir_ph = l_fir_c_d(1:end,:);
l_fir_amp = l_fir_c_d(1:end,:);
l_fir_amp_lnp = l_fir_amp;

for i = 1:size(l_sig,2)
    if strcmp(l_filter_flag ,'fir')
        % CIC Fir
        l_fir_c(:,i)      = filter(l_fir,l_mix_c(:,i)); % fir
        l_fir_s(:,i)      = filter(l_fir,l_mix_s(:,i)); % fir
%         l_fir_c_d(:,i)    = l_fir_c(l_wind_size:l_wind_size:end,i); % decimating
%         l_fir_s_d(:,i)    = l_fir_s(l_wind_size:l_wind_size:end,i); % decimating
        l_temp = filter(l_dec,l_fir_c(:,i)); % decimating
        l_fir_c_d(:,i)    = l_temp(1:size(l_fir_c_d,1));
        l_temp = filter(l_dec,l_fir_s(:,i)); % decimating
        l_fir_s_d(:,i)    = l_temp(1:size(l_fir_s_d,1));
        
    elseif strcmp(l_filter_flag ,'cic')
        % Fir and notch

        b_l = fir1(10, s_f_fofb/(2*f_adc));
        l_m = 1;          % Differential delays in the filter.
        l_n = 1;          % Filter sections
        l_r = l_wind_size;  % Decimation factor
        l_hm          = mfilt.cicdecim(l_r,l_m,l_n);             % Expects 16-bit input by default.

        l_fp_max = 2.5;   % max number in float point conversion
        l_fir_s_d_32  = filter(l_hm, int32(l_mix_s*(2^16)/l_fp_max));
        l_fir_c_d_32  = filter(l_hm, int32(l_mix_c*(2^16)/l_fp_max));
        l_fir_s_d     = double(l_fir_s_d_32)*l_fp_max/(2^16);  % converts back to float
        l_fir_c_d     = double(l_fir_c_d_32)*l_fp_max/(2^16);  % converts back to float
    end

    % ---
    % Cordic
    [l_fir_ph0(:,i), l_fir_amp0(:,i)] = cart2pol(l_fir_c_d(:,i), l_fir_s_d(:,i));

    l_fir_amp(:,i) = l_fir_amp0(1:end,i);   % correcting cic beginning
    l_fir_ph(:,i) = l_fir_ph0(1:end,i);     % correcting cic beginning
    
    % ---
    % Spike Removal
    if l_spk_rm
        l_amp_mean(i) = mean(l_fir_amp(:,1));
        l_fir_amp_temp = l_fir_amp(:,i);
        l_fir_amp_temp(l_fir_amp_temp >= l_amp_mean(i)*1.1) = l_amp_mean(i);
        l_fir_amp_temp(l_fir_amp_temp <= l_amp_mean(i)*0.9) = l_amp_mean(i);
        l_fir_amp(:,i) = l_fir_amp_temp;
    end
    
    % ---
    % Low-Pass filter
%     bl = fir1(14, 2000/(f_adc/l_wind_size/2));
    libera_lp_n = [0.0803696449793283 0.184926235006257 0.479755193683063 0.479755193683063 0.184926235006257 -0.0803696449793283]; %previously calculated with FDAtool
    % ---
    % Notch filter

    %[bn,an] = iirnotch((l_f_sw/4)/(f_adc/l_wind_size/2),(l_f_sw/4)/(f_adc/l_wind_size/2)/10);
    l_fir_amp_lpn(:,i) = filter(libera_lp_n,1,l_fir_amp(:,i));
end

    
    % Filtered Parameters
    l_f_fir   = round(f_adc/(l_wind_size));                    % updating decimated frequency
    l_t_fir   = linspace(0,max(l_t),length(l_fir_amp));

% ---
% Calculating X and Y - D/S and P D/S

[l_cic_xy_pds]        = calcpos_pds(l_fir_amp_lpn(30:end,:),Kx);       % PDS
[l_cic_xy_ds]         = calcpos_ds(l_fir_amp_lpn(30:end,:),Kx);        % DS

% picking "x" values
l_cic_x_pds = l_cic_xy_pds(:,1) - mean(l_cic_xy_pds(:,1));
l_cic_x_ds = l_cic_xy_ds(:,1) - mean(l_cic_xy_ds(:,1));

% Plotting

% ---
% Comparing integrated noise

% PSD calculation
l_wind    = rectwin(l_n_sw);
l_nover   = l_n_sw/2;
l_w       = l_n_sw;

l_cic_x_pds_psd       = pwelch(l_cic_x_pds,l_wind,l_nover,l_w,l_f_fir);

l_cic_x_ds_psd       = pwelch(l_cic_x_ds,l_wind,l_nover,l_w,l_f_fir);

% Integrating signals

l_f_step = l_f_fir/(2*length(l_cic_x_pds_psd)); % frequency/Hz step

l_cic_x_pds_psd_int       = cumtrapz(l_cic_x_pds_psd)*l_f_step;
l_cic_x_ds_psd_int       = cumtrapz(l_cic_x_ds_psd)*l_f_step;

%% Plotting results

s_fx_cic = linspace(0,s_f_cic/2,length(s_cic_x_pds_psd_int)); % frequency array
l_fx_cic = linspace(0,l_f_fir/2,length(l_cic_x_pds_psd_int)); % frequency array


figure
subplot(2,1,1)
loglog(s_fx_cic,sqrt(s_cic_x_pds_psd_int),'*-',s_fx_cic,sqrt(s_cic_x_ds_psd_int),'d-'); % integrated RMS
title('Squared Integrated Welch Power Spectral Density Estimate - Sirius');
ylabel('Integrated RMS noise [nm]')
xlabel('Frequency [Hz]')
legend('PDS - Arithmetic Mean','DS - Arithmetic Mean','location','best')
grid on

subplot(2,1,2)
loglog(l_fx_cic,sqrt(l_cic_x_pds_psd_int),'*-',l_fx_cic,sqrt(l_cic_x_ds_psd_int),'d-'); % integrated RMS
title('Squared Integrated Welch Power Spectral Density Estimate - Libera');
ylabel('Integrated RMS noise [nm]')
xlabel('Frequency [Hz]')
legend('PDS - Arithmetic Mean','DS - Arithmetic Mean','location','best')
grid on

figure
subplot(2,1,1)
loglog(s_fx_cic,sqrt(s_cic_x_pds_psd),'*-',s_fx_cic,sqrt(s_cic_x_ds_psd),'d-'); % integrated RMStitle('PSD Integrated - Noise')
title('Squared Welch Power Spectral Density Estimate - Sirius');
ylabel('RMS noise [nm]')
xlabel('Frequency [Hz]')
legend('PDS - Arithmetic Mean','DS - Arithmetic Mean','location','best')
grid on

subplot(2,1,2)
loglog(l_fx_cic,sqrt(l_cic_x_pds_psd),'*-',l_fx_cic,sqrt(l_cic_x_ds_psd),'d-'); % integrated RMStitle('PSD Integrated - Noise')
title('Squared Welch Power Spectral Density Estimate - Libera');
ylabel('RMS noise [nm]')
xlabel('Frequency [Hz]')
legend('PDS - Arithmetic Mean','DS - Arithmetic Mean','location','best')
grid on


% ---
% Finding x and y calculated values - Sirius

xy_ds_sirius = mean(s_cic_xy_ds);          % Sirius DS
xy_am_ds_sirius = mean(s_cic_xy_ds);       % Sirius DS am

xy_pds_sirius = mean(s_cic_xy_pds);        % Sirius PDS
xy_am_pds_sirius = mean(s_cic_xy_pds);     % Sirius PDS am

% sirius_string = ['Real Value       | ';'xy_ds_sirius     | '; 'xy_am_ds_sirius  | '; 'xy_pds_sirius    | '; 'xy_am_pds_sirius | '];
% sirius_xy     = [beam_xy; xy_ds_sirius ;xy_am_ds_sirius; xy_pds_sirius; xy_am_pds_sirius];
sirius_string = ['Real Value       | ';'xy_ds_sirius     | '; 'xy_pds_sirius    | '];
sirius_xy     = [beam_xy; xy_ds_sirius; xy_pds_sirius];

disp('Method           |   X (nm)            Y (nm)  ');
disp('-----------------------------------------------');
disp([sirius_string num2str(sirius_xy)]);

% ---
% Finding x and y calculated values - Libera

xy_ds_libera = mean(l_cic_xy_ds);             % Libera DS
xy_am_ds_libera = mean(l_cic_xy_ds);       % Libera DS am

xy_pds_libera = mean(l_cic_xy_pds);           % Libera PDS
xy_am_pds_libera = mean(l_cic_xy_pds);     % Libera PDS am

% libera_string = ['Real Value       | ';'xy_ds_libera     | '; 'xy_am_ds_libera  | '; 'xy_pds_libera    | '; 'xy_am_pds_libera | '];
% libera_xy     = [beam_xy; xy_ds_libera ;xy_am_ds_libera; xy_pds_libera; xy_am_pds_libera];
libera_string = ['Real Value       | '; 'xy_ds_libera     | '; 'xy_pds_libera    | '];
libera_xy     = [beam_xy; xy_ds_libera; xy_pds_libera];


%disp('Method           |   X (nm)            Y (nm)  ');
disp('-----------------------------------------------');
disp([libera_string num2str(libera_xy)]);