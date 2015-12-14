speech = wavread('speech.wav');
qout = zeros(length(speech),8);
mse = zeros(8,1);
for ii = 1:8
    [qout(:,ii),mse(ii)] = quantize(speech,ii);
end