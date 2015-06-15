fs = 8000;
y = zeros(8192,23);

figure(); hold on;
for k=1:23
    y(:,k) = load(sprintf('channel%d.txt',k));
    resp = abs(fft(y(:,k)));
    semilogx((0:(length(y)/2-1))*(fs/(2*length(y))),20*log10(resp(end/2+1:end)));
end

x = zeros(8192,1); x(1) = 1.0;
load cochlear_fb_8kHz_23_150_8192.mat

for k=1:23
    resp = abs(fft(G(:,k)));
    semilogx((0:(length(y)/2-1))*(fs/(2*length(y))),20*log10(resp(end/2+1:end)), 'r--');
end

hold off;