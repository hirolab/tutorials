%% Compare numerical derivatives in the frequency domain
addpath('..\..\mocap-utils\MATLAB')  % import deriv_sgolay

N = 10000;
T = 1 / 183;

x = randn(N,1);
rows = 10 : (N - 10);

xd_fw = (x(rows + 1) - x(rows)) / (1*T); 
xd_bw = (x(rows) - x(rows - 1)) / (1*T);
xd_ct = (x(rows + 1) - x(rows - 1)) / (2*T);
xd_sg = deriv_sgolay(x, 1/T, [5, 7]); xd_sg = xd_sg(rows);
x = x(rows);

[Txy, freq] = tfestimate(x, [xd_ct, xd_fw, xd_bw, xd_sg], hann(1000), 500, 1000, 1/T);

figure
h1 = subplot(221); semilogx(freq, 20*log10(abs([Txy, 1j*2*pi*freq])), 'linewidth', 1.5)
grid on, title('Magnitude'), ylabel('dB'), xlabel('Frequency [Hz]')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')
legend('forward', 'backward', 'central', 'Savitzky-Golay [5, 7] @ 183Hz', 'nominal', 'location', 'southeast')

h2 = subplot(223); semilogx(freq, angle([Txy, 1j*2*pi*freq])*180/pi, 'linewidth', 1.5)
grid on, title('Phase'), ylabel('deg'), xlabel('Frequency [Hz]')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')

h3 = subplot(222); semilogx(freq, 20*log10(abs(Txy./(1j*2*pi*freq))), 'linewidth', 1.5)
grid on, title('Magnitude Error'), ylabel('dB'), xlabel('Frequency [Hz]')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')

h4 = subplot(224); semilogx(freq, angle(Txy./(1j*2*pi*freq))*180/pi, 'linewidth', 1.5)
grid on, title('Phase Error'), ylabel('deg'), xlabel('Frequency [Hz]')
linkaxes([h1, h2, h3, h4], 'x'), xlim([1, 50]), 
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')

sgtitle('Types of Numerical Derivatives')
set(gcf, 'position', [680         344        1010         754])

print(gcf,'img/numerical_derivative.png','-dpng','-r300');