
clear all
close all
clc

%% Parameters Settings

% ------------------------------------------------------------
% Signal

f_rf    = 499.65794480373e6;
f_adc   = 203/864*f_rf;                 % sampling frequency
f_if    = 52/864*f_rf;                  % data frequency after sampling
f_rev   = f_rf/864;
f_fofb  = f_rev/6;
f_sw    = f_fofb;                       % switching frequency;
np_sw   = round(f_adc/f_sw);
n_sw    = 1000;                         % number of switching periods
np      = np_sw*n_sw;                   % number of points in data vector

% T       = [1.0 1.35 1.7 1.4];           % gains  T1, T2, T3 and T4 (Pairs are T1/T3, T2/T4)
T       = [1.0 1 1 1];           % gains  T1, T2, T3 and T4 (Pairs are T1/T3, T2/T4)
phases  = [-1.21 -0.70 -1.17 -1.00]*0;    % phase delays in rad
size_exp = 0;                         % size of the exponential, representing the switching transition

noise_flag = 0;                         % 1 for "noise", 0 for "without noise"

beam_xy = [10 10]*1e3;                  % position of the beam in  nm, [x y]
button_r = 3e6;                         % button radius in nm
chamber_r = 12e6;                       % chamber radius in nm
Kx = 8.5785e6;                          % converts to nm

% ------------------------------------------------------------
% Sirius

s_filter_flag = 'fir';                  % 'fir' or 'cic'

%% Creating signals - Sirius

% ---
% Creating time signal and frequency array using the sampling frequency

t       = (0:np-1)*(1/f_adc);
fx      = linspace(0,f_adc,np);

% Line gains specified on the beggining

% ---
% Data signal

abcd = pos2abcd(beam_xy,button_r,chamber_r); % obtaining the antennas intensity for the considered position

% phase delays specified in on the beginning

for i = 0:np/np_sw*2-1
    if rem(i,2) == 0 % even, delay of its channel
        s_sig(i*np_sw/2+1:((i+1)*np_sw)/2,1) = abcd(1) * sin(2*pi*f_rf*t(i*np_sw/2+1:((i+1)*np_sw)/2)+phases(1)); % Signal for line 1
        s_sig(i*np_sw/2+1:((i+1)*np_sw)/2,2) = abcd(2) * sin(2*pi*f_rf*t(i*np_sw/2+1:((i+1)*np_sw)/2)+phases(2)); % Signal for line 2
        s_sig(i*np_sw/2+1:((i+1)*np_sw)/2,3) = abcd(3) * sin(2*pi*f_rf*t(i*np_sw/2+1:((i+1)*np_sw)/2)+phases(3)); % Signal for line 3
        s_sig(i*np_sw/2+1:((i+1)*np_sw)/2,4) = abcd(4) * sin(2*pi*f_rf*t(i*np_sw/2+1:((i+1)*np_sw)/2)+phases(4)); % Signal for line 4
    else % odd, delay of opposing channel
        s_sig(i*np_sw/2+1:((i+1)*np_sw)/2,1) = abcd(1) * sin(2*pi*f_rf*t(i*np_sw/2+1:((i+1)*np_sw)/2)+phases(3)); % Signal for line 1
        s_sig(i*np_sw/2+1:((i+1)*np_sw)/2,2) = abcd(2) * sin(2*pi*f_rf*t(i*np_sw/2+1:((i+1)*np_sw)/2)+phases(4)); % Signal for line 2
        s_sig(i*np_sw/2+1:((i+1)*np_sw)/2,3) = abcd(3) * sin(2*pi*f_rf*t(i*np_sw/2+1:((i+1)*np_sw)/2)+phases(1)); % Signal for line 3
        s_sig(i*np_sw/2+1:((i+1)*np_sw)/2,4) = abcd(4) * sin(2*pi*f_rf*t(i*np_sw/2+1:((i+1)*np_sw)/2)+phases(2)); % Signal for line 4
    end
end

% ---
% Adding noise
if noise_flag
    noise   = randn(np,4)*1e-3;
    s_sig     = s_sig + noise;
end

% ---
% Creating Unswitched signal

s_sw0(:,1)   = square(2*pi*f_fofb*t(1:np_sw))*(T(1)-T(3))/2 + T(3) +(T(1)-T(3))/2;           % Signal for line 1
s_sw0(1:size_exp,1)   = (T(3)-T(1))*exp(-linspace(0,6,size_exp)) + T(1);                       % Transition 1
s_sw0(np_sw/2+1:np_sw/2+size_exp,1)   = (T(1)-T(3))*exp(-linspace(0,6,size_exp)) + T(3);   % Transition 2

s_sw0(:,2)   = square(2*pi*f_fofb*t(1:np_sw))*(T(2)-T(4))/2 + T(4) +(T(2)-T(4))/2;           % Signal for line 2
s_sw0(1:size_exp,2)   = (T(4)-T(2))*exp(-linspace(0,6,size_exp)) + T(2);                       % Transition 1
s_sw0(np_sw/2+1:np_sw/2+size_exp,2)   = (T(2)-T(4))*exp(-linspace(0,6,size_exp)) + T(4);   % Transition 2

s_sw0(:,3)   = square(2*pi*f_fofb*t(1:np_sw)+pi)*(T(1)-T(3))/2 + T(3) +(T(1)-T(3))/2;        % Signal for line 3
s_sw0(1:size_exp,3)   = (T(1)-T(3))*exp(-linspace(0,6,size_exp)) + T(3);                       % Transition 1
s_sw0(np_sw/2+1:np_sw/2+size_exp,3)   = (T(3)-T(1))*exp(-linspace(0,6,size_exp)) + T(1);   % Transition 2

s_sw0(:,4)   = square(2*pi*f_fofb*t(1:np_sw)+pi)*(T(2)-T(4))/2 + T(4) +(T(2)-T(4))/2;        % Signal for line 4
s_sw0(1:size_exp,4)   = (T(2)-T(4))*exp(-linspace(0,6,size_exp)) + T(4);                       % Transition 1
s_sw0(np_sw/2+1:np_sw/2+size_exp,4)   = (T(4)-T(2))*exp(-linspace(0,6,size_exp)) + T(2);   % Transition 2

s_sw = repmat(s_sw0,n_sw,1);

s_sig_sw = zeros(np_sw*n_sw,4); % allocating variable
for i = 1:size(s_sig,2)
    s_sig_sw(:,i)  = s_sig(:,i) .* s_sw(:,i);
end

% ---
% Plotting switched signal
figure
plot(1:np_sw,s_sw0(:,1),1:np_sw,square(2*pi*f_fofb*t(1:np_sw))*(T(1)-T(3))/2 + T(3) +(T(1)-T(3))/2,'r');
title('Switched Signal - Channel 1 (A and C) Amplitudes');
ylabel('Amplitude')
xlabel('Time')
grid on
legend('Attenuated','Ideal','location','best')
axis([1 size(s_sw0,1) min(s_sw0(:,1))*0.95 max(s_sw0(:,1))*1.05]);

% ---
% Plotting switched signal
figure
plot(t,s_sig_sw(:,1),t,s_sw(:,1),'r');
title('Unswitched Signal - Channel 1 (A and C)');
ylabel('Amplitude')
xlabel('Time [s]')
legend('Unswitched Signal','Switching','location','best')
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
    s_cos       = cos(2*pi*f_fofb*t);
    s_sin       = sin(2*pi*f_fofb*t);
    s_mix_c(:,i)  = s_sig_sw(:,i).*s_cos';
    s_mix_s(:,i)  = -1*s_sig_sw(:,i).*s_sin';

    % ---
    % Filter
    s_wind_size   = round(f_adc/(f_fofb));      % size of the windows

    if strcmp(s_filter_flag ,'fir')
        % Fir
        s_cic             = dfilt.dffir(((1/(s_wind_size)*ones(1,s_wind_size)))); % creating fir equivalent to cic
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

    % ---
    % Cordic
    [s_cic_ph0(:,i), s_cic_amp0(:,i)] = cart2pol(s_cic_c_d(:,i), s_cic_s_d(:,i));

    s_cic_amp(:,i) = s_cic_amp0(2:end,i); % correcting cic beginning

    % ---
    % Mean values
    s_cic_amp_am(:,i) = (s_cic_amp(:,i));

    % Filtered Parameters
    s_f_cic   = round(f_adc/(s_wind_size));                    % updating decimated frequency
    s_t_cic   = linspace(0,max(t),length(s_cic_amp));

end

% ---
% Calculating X and Y - D/S and P D/S

[s_cic_xy_pds]        = calcpos_pds(s_cic_amp(1:end,:),Kx);       % PDS - na
[s_cic_xy_ds]         = calcpos_ds(s_cic_amp(1:end,:),Kx);        % DS  - na

[s_cic_xy_am_pds]     = calcpos_pds(s_cic_amp_am(1:end,:),Kx);    % PDS - am
[s_cic_xy_am_ds]      = calcpos_ds(s_cic_amp_am(1:end,:),Kx);     % DS  - am

% picking "x" values
s_cic_x_pds_ac = s_cic_xy_pds(:,1) - mean(s_cic_xy_pds(:,1));
s_cic_x_ds_ac = s_cic_xy_ds(:,1) - mean(s_cic_xy_ds(:,1));

s_cic_x_am_pds_ac = s_cic_xy_am_pds(:,1) - mean(s_cic_xy_am_pds(:,1));
s_cic_x_am_ds_ac = s_cic_xy_am_ds(:,1) - mean(s_cic_xy_am_ds(:,1));

% Plotting

% ---
% Comparing integrated noise

% PSD calculation
s_wind    = rectwin(n_sw/2);
s_nover   = n_sw/4;
s_w       = n_sw/2;

s_cic_x_pds_psd       = pwelch(s_cic_x_pds_ac,s_wind,s_nover,s_w,s_f_cic);
s_cic_x_am_pds_psd    = pwelch(s_cic_x_am_pds_ac,s_wind,s_nover,s_w,s_f_cic);

s_cic_x_ds_psd       = pwelch(s_cic_x_pds_ac,s_wind,s_nover,s_w,s_f_cic);
s_cic_x_am_ds_psd    = pwelch(s_cic_x_am_pds_ac,s_wind,s_nover,s_w,s_f_cic);

% Integrating signals

s_f_step = s_f_cic/(2*length(s_cic_x_pds_psd)); % frequency/Hz step

s_cic_x_pds_psd_int       = cumtrapz(s_cic_x_pds_psd)*s_f_step;
s_cic_x_am_pds_psd_int    = cumtrapz(s_cic_x_am_pds_psd)*s_f_step;

s_cic_x_ds_psd_int       = cumtrapz(s_cic_x_ds_psd)*s_f_step;
s_cic_x_am_ds_psd_int    = cumtrapz(s_cic_x_am_ds_psd)*s_f_step;

% % Erasing DC level
% cic_x_pds_psd_int       = cic_x_pds_psd_int - cic_x_pds_psd_int(1);
% cic_x_am_pds_psd_int    = cic_x_am_pds_psd_int - cic_x_am_pds_psd_int(1);
% cic_x_gm_pds_psd_int    = cic_x_gm_pds_psd_int - cic_x_gm_pds_psd_int(1);
%
% cic_x_ds_psd_int       = cic_x_ds_psd_int - cic_x_ds_psd_int(1);
% cic_x_am_ds_psd_int    = cic_x_am_ds_psd_int - cic_x_am_ds_psd_int(1);
% cic_x_gm_ds_psd_int    = cic_x_gm_ds_psd_int - cic_x_gm_ds_psd_int(1);


%% Plotting results

s_fx_cic = linspace(0,s_f_cic/2,length(s_cic_x_pds_psd_int)); % frequency array

figure
loglog(s_fx_cic,sqrt(s_cic_x_pds_psd_int),'*-',s_fx_cic,sqrt(s_cic_x_ds_psd_int),'d-'); % integrated RMS
title('Squared Integrated Welch Power Spectral Density Estimate - Sirius');
ylabel('Integrated RMS noise [nm]')
xlabel('Frequency [Hz]')
legend('PDS - Arithmetic Mean','DS - Arithmetic Mean','location','best')
grid on

figure
loglog(s_fx_cic,sqrt(s_cic_x_pds_psd),'*-',s_fx_cic,sqrt(s_cic_x_ds_psd),'d-'); % integrated RMStitle('PSD Integrated - Noise')
title('Squared Welch Power Spectral Density Estimate - Sirius');
ylabel('RMS noise [nm]')
xlabel('Frequency [Hz]')
legend('PDS - Arithmetic Mean','DS - Arithmetic Mean','location','best')
grid on

% ---
% Finding x and y calculated values - Sirius

xy_ds_sirius = mean(s_cic_xy_ds);             % Sirius DS
xy_am_ds_sirius = mean(s_cic_xy_am_ds);       % Sirius DS am

xy_pds_sirius = mean(s_cic_xy_pds);           % Sirius PDS
xy_am_pds_sirius = mean(s_cic_xy_am_pds);     % Sirius PDS am

% sirius_string = ['Real Value       | ';'xy_ds_sirius     | '; 'xy_am_ds_sirius  | '; 'xy_pds_sirius    | '; 'xy_am_pds_sirius | '];
% sirius_xy     = [beam_xy; xy_ds_sirius ;xy_am_ds_sirius; xy_pds_sirius; xy_am_pds_sirius];
sirius_string = ['Real Value       | ';'xy_ds_sirius     | '; 'xy_pds_sirius    | '];
sirius_xy     = [beam_xy; xy_ds_sirius; xy_pds_sirius];

disp('Method           |   X (nm)            Y (nm)  ');
disp('-----------------------------------------------');
disp([sirius_string num2str(sirius_xy)]);

%%

alfa1 = T(1)+T(3);
alfa2 = T(2)+T(4);

a = abcd(1);
b = abcd(2);
c = abcd(3);
d = abcd(4);