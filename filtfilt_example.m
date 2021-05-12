fs = 2000; %rate that data is being recorded at
t = 1/fs:1/fs:1;

n1 = sin(2*pi*10*t); %vibrations from skin signal
s1 = cos(2*pi*35*t);  %EMG data
n2 = randn(length(t),1)';  %noise not good

s = n1+s1+n2;
s_ideal = s1; 

figure;plot(t,s,t,s_ideal)

fn = 2/fs;

[b a] = fir1(100,[30*fn 50*fn])

freqz(b,1,512)


y = filtfilt(b,a,s);
figure;plot(t,s,t,s_ideal,t,y,'LineWidth',2)
legend('noisy','ideal','filtered')



