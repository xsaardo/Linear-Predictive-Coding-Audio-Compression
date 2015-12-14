[speech, Fs] = wavread('asknot.wav');

% Parameters to Set
alpha = 3;
r = 3;

L = 2^r;
m = mean(speech);
sigma = std(speech);

speech(speech >= sigma*alpha) = sigma*alpha;

q = (max(speech) - min(speech))/L;

yq = round(speech/q)*q;
