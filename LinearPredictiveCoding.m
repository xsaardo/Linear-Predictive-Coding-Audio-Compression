clc;
close all;
clear all;

%% User Parameters(Set as necessary)
alpha = 3;          % same alpha for both speech truncation and residual quantization
BIT = 8;            % Quantization Rate
resFlag = true;     % Flag for quantizing residuals
manResFlag = false; % Flag for manual residual calculation
paramFlag = true;   % Flag for coefficient quantization
alpha_coeff = 3;    % alpha for coefficient quantization
BLOCK_LENGTH = 160; % Block length
ORDER = 10;         % AR filter order

speech = audioread('futuresound.wav');              % Read in speech file
qout = zeros(length(speech),8);                     % Initialize output file
mse = zeros(8,1);                                   % Initialize MSE vector for direct quantization

%% Direct Quantization and Plots
figure;
for a = 1:8             % For alpha ranging from 1 to 8
    for ii = 1:8        % For quantization rate ranging from 1 to 8
        [qout(:,ii),mse(ii)] = quantize(speech,ii,a);   % Quantize speech and save MSE in mse vector
    end
    hold on
    plot(mse);
end
xlabel('Quantization Rate r');
ylabel('MSE');
legend('a=1','a=2','a=3','a=4','a=5','a=6');

%% Outlier Truncation Preprocessing
m = mean(speech);
sigma = std(speech);
speech(speech >= m + sigma*alpha) = sigma*alpha + m;
speech(speech <= m - sigma*alpha) = -sigma*alpha + m;

%% LPC Parameters
numblocks = ceil(length(speech)/BLOCK_LENGTH);
ytotal = [];                                % initialize output sound matrix
mse_quantres = [];                          % initialize mse matrix

for coeff_bit = [4 8 9]         % For coefficient quantization rates 4,8,9
    for res_bit = 1:8           % For residual quantization rates 1 to 8
        
        % Cutting off last block 
        speech = speech(1:(numblocks-1)*BLOCK_LENGTH);
        
        % Residual/Coefficient Calculations
        mse_rec = zeros(numblocks,1);
        for ii = 0:(numblocks-2)
            
            % First Block Case
            if ii == 0
                A = toeplitz([0;speech(1:BLOCK_LENGTH-1)],zeros(ORDER,1));
                param = A\speech(1:BLOCK_LENGTH); % Parameter Calculations
                residual1 = speech(1:BLOCK_LENGTH) - A*param; % Residual Calculations
                if paramFlag
                    param = quantize(param,coeff_bit,alpha_coeff);
                end
                if resFlag
                    residual1 = quantize(residual1,res_bit,alpha);
                end
                
                % Manual Residual Calculation
                if manResFlag
                    residual2 = zeros(BLOCK_LENGTH,1);
                    for jj = 1:BLOCK_LENGTH
                        if jj == 1
                            residual2(jj) = speech(jj);
                        elseif jj <= ORDER && jj > 1
                            residual2(jj) = speech(jj) - param'*[flipud(speech(1:jj-1));zeros(ORDER - jj + 1,1)];
                        else
                            residual2(jj) = speech(jj) - param'*flipud(speech(jj-ORDER:jj-1));
                        end
                    end
                    
                    % MSE Between Residuals
                    mse_res(ii+1) = (residual1-residual2)'*(residual1-residual2)/BLOCK_LENGTH;
                end
                
                
                % Reconstruction
                y = zeros(BLOCK_LENGTH,1);
                for kk = 1:BLOCK_LENGTH
                    if kk <= ORDER
                        y(kk) = param'*[flipud(y(1:kk-1));zeros(ORDER-kk+1,1)] + residual1(kk);
                    else
                        y(kk) = param'*flipud(y(kk-ORDER:kk-1)) + residual1(kk);
                    end
                end
                
                % Reconstruction MSE
                mse_rec(ii+1) = (y-speech(1:BLOCK_LENGTH))'*(y-speech(1:BLOCK_LENGTH))/BLOCK_LENGTH;
                
                % Rest of the Blocks
            else
                prevSeg = speech((ii-1)*BLOCK_LENGTH + 1:ii*BLOCK_LENGTH);
                speechSeg = speech(ii*BLOCK_LENGTH + 1:(ii + 1)*BLOCK_LENGTH);
                A = toeplitz([prevSeg(end);speechSeg(1:BLOCK_LENGTH-1)],flipud(prevSeg(end-ORDER+1:end)));
                param = A\speechSeg;
                residual1 = speechSeg - A*param;
                if paramFlag
                    param = quantize(param,coeff_bit,alpha_coeff);
                end
                if resFlag
                    residual1 = quantize(residual1,res_bit,alpha);
                end
                
                % Manual Residual Calculation
                if manResFlag
                    for jj = 1:BLOCK_LENGTH
                        if jj == 1
                            residual2(jj) = speechSeg(jj) - param'*flipud(prevSeg(end-ORDER+1:end));
                        elseif jj <= ORDER && jj > 1
                            residual2(jj) = speechSeg(jj) - param'*[flipud(speechSeg(1:jj-1));flipud(prevSeg(end-ORDER + jj:end))];
                        else
                            residual2(jj) = speechSeg(jj) - param'*flipud(speechSeg(jj-ORDER:jj-1));
                        end
                    end
                    % MSE Between Residuals
                    mse_res(ii+1) = (residual1-residual2)'*(residual1-residual2)/BLOCK_LENGTH;
                end
                
                
                % Reconstruction
                yblock = zeros(BLOCK_LENGTH,1);
                for jj = 1:BLOCK_LENGTH
                    if jj <= ORDER
                        yblock(jj) = param'*[flipud(yblock(1:jj-1));flipud(y(end-ORDER+jj:end))] + residual1(jj);
                    else
                        yblock(jj) = param'*flipud(yblock(jj-ORDER:jj-1)) + residual1(jj);
                    end
                end
                y = [y;yblock];
                
                % Reconstruction MSE
                mse_rec(ii+1) = (yblock-speechSeg)'*(yblock-speechSeg)/BLOCK_LENGTH;
                
            end
        end
        
        % Matrix of MSE values (columns refer to coefficient quantization rates/rows
        % correspond to residual quantization rates)
        mse_quantres(res_bit,coeff_bit) = (speech-y)'*(speech-y)/length(speech); 
        
        % Output speech matrix (First 8 columns refer to coefficient
        % quantization rate 4, second 8 rate 8, last 8 rate 9)
        % Each set of 8 columns refers to 8 residual quantization rates
        ytotal = [ytotal y];
    end
end