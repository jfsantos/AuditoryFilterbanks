fs = 16000;
%y = zeros(512,23);
figure(); hold on;
for k=1:23
    %y(:,k) = load(sprintf('channel%d.txt',k));
    resp = abs(fft(y(:,k)));
    semilogx((0:(length(y)/2-1))*(fs/(2*length(y))),20*log10(resp(end/2+1:end)));
end

hold off;