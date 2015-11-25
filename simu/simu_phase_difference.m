
clear all
close all
clc

%% Creating signals

% ---
% System parameters

f_rf    = 499.65794480373e6;    % sirius data frequency
f_adc   = 203/864*f_rf;         % sampling frequency
f_if    = 52/864*f_rf;          % data frequency after sampling

f_rev   = f_rf/864;
f_fofb  = f_rev/10;              % switching frequency;
np_fofb = f_adc/f_fofb;
np      = np_fofb*1000;  % number of points in data vector

% ---
% Creating switching signal and frequency array

t       = (0:np-1)*(1/f_adc);
fx      = linspace(0,f_adc,np);

T       = [3679 3680 3603 3610/2]; % gains  T1, T2, T3 and T4

sw(:,1)   = square(2*pi*f_fofb*t)*(T(1)-T(3))/2 + T(3) +(T(1)-T(3))/2;        % Signal for line 1
sw(:,3)   = square(2*pi*f_fofb*t + pi)*(T(1)-T(3))/2 + T(3) +(T(1)-T(3))/2;   % Signal for line 3
sw(:,2)   = square(2*pi*f_fofb*t)*(T(2)-T(4))/2 + T(4) +(T(2)-T(4))/2;        % Signal for line 2
sw(:,4)   = square(2*pi*f_fofb*t + pi)*(T(2)-T(4))/2 + T(4) +(T(2)-T(4))/2;   % Signal for line 4

% ---
% Data signal

phases  = [pi/60 0 pi/2 pi/90]; % phase delays in rad

for i = 0:np/np_fofb*2-1
    if rem(i,2) % odd, no delay
        sig(i*np_fofb/2+1:((i+1)*np_fofb)/2,1) = sin(2*pi*f_rf*t(i*np_fofb/2+1:((i+1)*np_fofb)/2)); % Signal for line 1
        sig(i*np_fofb/2+1:((i+1)*np_fofb)/2,2) = sin(2*pi*f_rf*t(i*np_fofb/2+1:((i+1)*np_fofb)/2)); % Signal for line 2
        sig(i*np_fofb/2+1:((i+1)*np_fofb)/2,3) = sin(2*pi*f_rf*t(i*np_fofb/2+1:((i+1)*np_fofb)/2)); % Signal for line 3
        sig(i*np_fofb/2+1:((i+1)*np_fofb)/2,4) = sin(2*pi*f_rf*t(i*np_fofb/2+1:((i+1)*np_fofb)/2)); % Signal for line 4
    else % even, delayed
        sig(i*np_fofb/2+1:((i+1)*np_fofb)/2,1) = sin(2*pi*f_rf*t(i*np_fofb/2+1:((i+1)*np_fofb)/2)+phases(1)); % Signal for line 1
        sig(i*np_fofb/2+1:((i+1)*np_fofb)/2,2) = sin(2*pi*f_rf*t(i*np_fofb/2+1:((i+1)*np_fofb)/2)+phases(2)); % Signal for line 2
        sig(i*np_fofb/2+1:((i+1)*np_fofb)/2,3) = sin(2*pi*f_rf*t(i*np_fofb/2+1:((i+1)*np_fofb)/2)+phases(3)); % Signal for line 3
        sig(i*np_fofb/2+1:((i+1)*np_fofb)/2,4) = sin(2*pi*f_rf*t(i*np_fofb/2+1:((i+1)*np_fofb)/2)+phases(4)); % Signal for line 4
    end
end

 carrier = sin(2*pi*f_rf*t);

% sig     = [carrier' carrier' carrier' carrier'];
noise   = randn(np,4)*1e-3;

% sig     = sig + noise;

% only for tests
sig_test0 = sig(:,1) .* sw(:,1);
sig_test1 = carrier' .* sw(:,1);

% ---
% Switching signals

for i = 1:size(sig,2)    
    sig_sw(:,i)  = sig(:,i) .* sw(:,i);
end

% Plotting switched signal
figure
plot(t,sig_sw(:,1),t,sw(:,1),'r');
title('Switched Signal with Noise - Channel 1');
ylabel('Amplitude')
xlabel('Time [s]')
legend('Switched Signal with noise','Switching','location','best')
grid on

%% Processing - Sirius

for i = 1:size(sig,2)
    
    % ---
    % Mixing signal
    s_cos       = cos(2*pi*f_fofb*t);
    s_sin       = sin(2*pi*f_fofb*t);
    mix_c(:,i)  = sig_sw(:,i).*s_cos';
    mix_s(:,i)  = -1*sig_sw(:,i).*s_sin';
    
    % ---
    % Filter
    wind_size   = round(f_adc/(2*f_fofb));      % size of the windows

    
%     % Fir   
%     cic             = dfilt.dffir(((1/(wind_size)*ones(1,wind_size)))); % creating fir equivalent to cic 
%     cic_c(:,i)      = filter(cic,mix_c(:,i)); % cic
%     cic_s(:,i)      = filter(cic,mix_s(:,i)); % cic
%     cic_c_d(:,i)    = cic_c(wind_size:wind_size:end,i); % decimating
%     cic_s_d(:,i)    = cic_s(wind_size:wind_size:end,i); % decimating
    
    % CIC
    m = 1;          % Differential delays in the filter.
    n = 1;          % Filter sections
    r = wind_size;  % Decimation factor
    hm          = mfilt.cicdecim(r,m,n);             % Expects 16-bit input by default.
    cic_s_d_32  = filter(hm, int32(mix_s));
    cic_c_d_32  = filter(hm, int32(mix_c));
    cic_s_d     = double(cic_s_d_32);
    cic_c_d     = double(cic_c_d_32);

    
    % ---
    % Cordic
    [cic_ph0(:,i), cic_amp0(:,i)] = cart2pol(cic_c_d(:,i), cic_s_d(:,i));
    
    cic_amp(:,i) = cic_amp0(2:end,i); % correcting cic beginning
    
    % ---
    % Mean values
    cic_amp_am(:,i) = amean(cic_amp(:,i));
    cic_amp_gm(:,i) = gmean(cic_amp(:,i));
    
    % Filtered Parameters
    f_cic   = round(f_adc/(wind_size));                    % updating decimated frequency
    t_cic   = linspace(0,max(t),length(cic_amp));
    
end

% ---
% Calculating X and Y - D/S and P D/S

Kx = 8.5785e6; % converts to nm

[cic_xy_pds]        = calcpos_pds(cic_amp(1:end,:),Kx);       % PDS - na
[cic_xy_ds]         = calcpos_ds(cic_amp(1:end,:),Kx);        % DS  - na

[cic_xy_am_pds]     = calcpos_pds(cic_amp_am(1:end,:),Kx);    % PDS - am
[cic_xy_am_ds]      = calcpos_ds(cic_amp_am(1:end,:),Kx);     % DS  - am

[cic_xy_gm_pds]     = calcpos_pds(cic_amp_gm(1:end,:),Kx);    % PDS - gm
[cic_xy_gm_ds]      = calcpos_ds(cic_amp_gm(1:end,:),Kx);     % DS  - gm

% picking "x" values
cic_x_pds = cic_xy_pds(:,1) - mean(cic_xy_pds(:,1));
cic_x_ds = cic_xy_ds(:,1) - mean(cic_xy_ds(:,1));

cic_x_am_pds = cic_xy_am_pds(:,1) - mean(cic_xy_am_pds(:,1));
cic_x_am_ds = cic_xy_am_ds(:,1) - mean(cic_xy_am_ds(:,1));

cic_x_gm_pds = cic_xy_gm_pds(:,1) - mean(cic_xy_gm_pds(:,1));
cic_x_gm_ds = cic_xy_gm_ds(:,1) - mean(cic_xy_gm_ds(:,1));

% Plotting 

% ---
% Comparing integrated noise

% PSD calculation
wind    = rectwin(1000);
nover   = 500;
w       = 1000;

cic_x_pds_psd       = pwelch(cic_x_pds,wind,nover,w,f_cic);
cic_x_am_pds_psd    = pwelch(cic_x_am_pds,wind,nover,w,f_cic);
cic_x_gm_pds_psd    = pwelch(cic_x_gm_pds,wind,nover,w,f_cic);

cic_x_ds_psd       = pwelch(cic_x_pds,wind,nover,w,f_cic);
cic_x_am_ds_psd    = pwelch(cic_x_am_pds,wind,nover,w,f_cic);
cic_x_gm_ds_psd    = pwelch(cic_x_gm_pds,wind,nover,w,f_cic);

% Integrating signals

f_step = f_cic/(2*length(cic_x_pds_psd)); % frequency/Hz step

cic_x_pds_psd_int       = cumtrapz(cic_x_pds_psd)*f_step;
cic_x_am_pds_psd_int    = cumtrapz(cic_x_am_pds_psd)*f_step;
cic_x_gm_pds_psd_int    = cumtrapz(cic_x_gm_pds_psd)*f_step;

cic_x_ds_psd_int       = cumtrapz(cic_x_ds_psd)*f_step;
cic_x_am_ds_psd_int    = cumtrapz(cic_x_am_ds_psd)*f_step;
cic_x_gm_ds_psd_int    = cumtrapz(cic_x_gm_ds_psd)*f_step;

% % Erasing DC level
% cic_x_pds_psd_int       = cic_x_pds_psd_int - cic_x_pds_psd_int(1);
% cic_x_am_pds_psd_int    = cic_x_am_pds_psd_int - cic_x_am_pds_psd_int(1);
% cic_x_gm_pds_psd_int    = cic_x_gm_pds_psd_int - cic_x_gm_pds_psd_int(1);
% 
% cic_x_ds_psd_int       = cic_x_ds_psd_int - cic_x_ds_psd_int(1);
% cic_x_am_ds_psd_int    = cic_x_am_ds_psd_int - cic_x_am_ds_psd_int(1);
% cic_x_gm_ds_psd_int    = cic_x_gm_ds_psd_int - cic_x_gm_ds_psd_int(1);

% Plotting results

fx_cic = linspace(0,f_cic/2,length(cic_x_pds_psd_int)); % frequency array

figure
loglog(fx_cic,sqrt(cic_x_pds_psd_int),'*-',fx_cic,sqrt(cic_x_am_pds_psd_int),'*-',fx_cic,sqrt(cic_x_gm_pds_psd_int),'*-',fx_cic,sqrt(cic_x_ds_psd_int),'d-',fx_cic,sqrt(cic_x_am_ds_psd_int),'d-',fx_cic,sqrt(cic_x_gm_ds_psd_int),'d-'); % integrated RMS
title('Squared Integrated Welch Power Spectral Density Estimate');
ylabel('Integrated RMS noise [nm]')
xlabel('Frequency [Hz]')
legend('PDS - No Averaged','PDS - Arithmetic Mean','PDS - Geometric Mean','DS - No Averaged','DS - Arithmetic Mean','DS - Geometric Mean','location','best')
grid on

figure
loglog(fx_cic,sqrt(cic_x_pds_psd),'*-',fx_cic,sqrt(cic_x_am_pds_psd),'*-',fx_cic,sqrt(cic_x_gm_pds_psd),'*-',fx_cic,sqrt(cic_x_ds_psd),'d-',fx_cic,sqrt(cic_x_am_ds_psd),'d-',fx_cic,sqrt(cic_x_gm_ds_psd),'d-'); % integrated RMStitle('PSD Integrated - Noise')
title('Squared Welch Power Spectral Density Estimate');
ylabel('RMS noise [nm]')
xlabel('Frequency [Hz]')
legend('PDS - No Averaged','PDS - Arithmetic Mean','PDS - Geometric Mean','DS - No Averaged','DS - Arithmetic Mean','DS - Geometric Mean','location','best')
grid on
