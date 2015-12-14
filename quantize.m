% Helper function for LPC that performs outlier truncation and quantization to r bits or 2^r distinct
% levels and calculates MSE as well

function [qout,mse] = quantize(input,r,alpha)
    L = 2^r;
    
    %Outlier Truncation
    m = mean(input);
    sigma = std(input);
    input(input >= m + sigma*alpha) = sigma*alpha + m;
    input(input <= m - sigma*alpha) = -sigma*alpha + m;
     
    % Quantization
    q = (max(input) - min(input))/L;
    qout = round(input/q)*q;
    
    % MSE Calculation
    mse = ((input - qout)'*(input - qout))/length(input);
end