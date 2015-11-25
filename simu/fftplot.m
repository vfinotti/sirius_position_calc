function fftplot(x,fs)

subplot(2,1,1)
plot(x)
title('Time Plot')
grid on

subplot(2,1,2)
plot(linspace(0,fs,length(x)),abs(fft(x)))
title('Frequency Plot')
grid on

end