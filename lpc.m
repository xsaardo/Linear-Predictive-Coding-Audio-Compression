% Project for ECE 174 (Linear/Nonlinear optimization) in which we compress speech by quantizing its
% residuals.  The residuals are calculated by subtracting the sample speech
% from an autoregressive moving average model.  

speech = wavread('asknot.wav');
BLOCK_LENGTH = 160;
ORDER = 10;
BIT = 1;
numblocks = ceil(length(speech)/BLOCK_LENGTH);

% Zero Padding Speech Signal
speech = [speech;zeros(BLOCK_LENGTH*numblocks - length(speech),1)];

%% Parameter Estimation
y = zeros(BLOCK_LENGTH,1);
for ii = 0:numblocks-1 
    if ii == 0
        A = toeplitz([0;speech(1:BLOCK_LENGTH-1)],zeros(ORDER,1));
        param = A\speech(1:BLOCK_LENGTH);
        residual = speech(1:BLOCK_LENGTH) - A*param;
        [residual,~] = quantize(residual,BIT,1);
        for kk = 1:BLOCK_LENGTH
            if kk <= ORDER
                y(kk) = param'*[flipud(y(1:kk-1));zeros(ORDER-kk+1,1)] + residual(kk);
            else
                y(kk) = param'*flipud(y(kk-ORDER:kk-1)) + residual(kk);
            end
        end
    else 
        prevSeg = speech((ii-1)*BLOCK_LENGTH + 1:ii*BLOCK_LENGTH);
        speechSeg = speech(ii*BLOCK_LENGTH + 1:(ii + 1)*BLOCK_LENGTH);
        A = toeplitz([prevSeg(end);speechSeg(1:BLOCK_LENGTH-1)],flipud(prevSeg(end-ORDER+1:end)));
        param = A\speechSeg;
        residual = speechSeg - A*param;
        [residual,~] = quantize(residual,BIT,1);
        
        % Manual Residual Calculation + MSE
        % error = speechSeg - (param'*toeplitz(flipud(prevSeg(end-ORDER+1:end)),[prevSeg(end);speechSeg(1:BLOCK_LENGTH-1)]))';
        % mean((residual - error).^2);
        
        % Reconstruction
        yblock = zeros(BLOCK_LENGTH,1);
        for jj = 1:BLOCK_LENGTH
            if jj <= ORDER
                yblock(jj) = param'*[flipud(yblock(1:jj-1));flipud(y(end-ORDER+jj:end))] + residual(jj);
            else
                yblock(jj) = param'*flipud(yblock(jj-ORDER:jj-1)) + residual(jj);
            end
        end
        y = [y;yblock];
    end
end