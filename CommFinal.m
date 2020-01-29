% Gian Angelo Tria
% ECE-408 | Wireless Communications
% 01/29/19 

% Goal: Design and simulate an end to end communication link with Encoding, Modulation, and Equalization in MATLAB.

% Objectives: 
% Part 1: Match AWGN performance for 4 and 16 QAM with no ISI. 
% Achieve BER of 10^-4 at 12 dB SNR using BPSK + Adaptive Equalizer over moderate ISI channel.
% Part 2: Achieve BER of 10^-6 at 12 dB SNR over moderate ISI channel using whatever means possible.
%% PART 1
% A skeleton BER script for a wireless link simulation
clear all;close all;clc

numIter = 1;  % The number of iterations of the simulation
nSym = 1000;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

M = [2 4 16];        % The M-ary number, 2 corresponds to binary modulation

% chan = 1;          % No channel
chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI
%chan = [0.227 0.460 0.688 0.460 0.227]';   % Not so invertible, severe ISI


% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);
berVec2 = zeros(numIter, lenSNR);

% Equalizer Object
% eqobj = lineareq(6, lms(0.01)); % Linear LMS
eqobj = lineareq(6, rls(1)); % Attempt at Linear RLS
% eqobj = dfe(6,7,lms(0.01)); %Attempt at Decision Feedback

for m = M
% Run the simulation numIter amount of times
for i = 1:numIter

    %bits = randi([0 1], 1, nSym*m);     % Generate random bits
    bits = randi(2,[nSym*log2(m), 1])-1; % Generate random bits 2
    % New bits must be generated at every iteration
    
   % If you increase the M-ary number, as you most likely will, you'll need to
    % convert the bits to integers. See the BIN2DE function
    % For binary, our MSG signal is simply the bits
      msg = convertB2M(bits, m);

    for j = 1:lenSNR % one iteration of the simulation at each SNR Value


        tx = qammod(msg,m);  % BPSK modulate the signal

        if isequal(chan,1)
            txChan = tx;
        else
            txChan = filter(chan,1,tx);  % Apply the channel.
        end
% First Attempt
%             weight = zeros(lenSNR,1);
%             x = zeros(lenSNR,1);
%             error = [];
%             for i = lenSNR+1:nSym*M
%                 x = txNoisy(i:-1:i-(lenSNR-1));
%                 y = x*weight(:,1);
%                 error(i) = bits(i) - y;
%                 weight = weight + u*x*error(i);
%             end

        % txNoisy = awgn(txChan,SNR_Vec(j)); % Add AWGN
        
        %Adding AWGN - Properly scaling noise
        if isequal(m,2)
            noise = 0;
        else
            noise =  round(10*log10(log2(m)));
        end
        txNoisy = awgn(txChan,SNR_Vec(j)+noise);
        txNoisy2 = equalize(eqobj,txNoisy,tx(1:100));
        
        rx = qamdemod(txNoisy,m); % Demodulate
        rx2 = qamdemod(txNoisy2,m); % Demodulate
        % Again, if M was a larger number, I'd need to convert my symbols
        % back to bits here.
        %rxMSG = rx;
        %rxMSG2 = rx2;
        rxMSG = convertM2B(rx, m)';
        rxMSG2 = convertM2B(rx2, m)';
        % Compute and store the BER for this iteration

        [~, berVec(i,j)] = biterr(bits, rxMSG);  % We're interested in the BER, which is the 2nd output of BITERR
        [~, berVec2(i,j)] = biterr(bits, rxMSG2); 
    end  % End SNR iteration
end      % End numIter iteration


% Compute and plot the mean BER
ber = mean(berVec,1);
ber2 = mean(berVec2,1);
figure
semilogy(SNR_Vec, ber) %Use Log scale for SNR

% Compute the theoretical BER for this scenario
berTheory = berawgn(SNR_Vec,'psk',2,'nondiff');
hold on
xlabel('E_b/N_0(dB)');  ylabel('BER');
semilogy(SNR_Vec,berTheory,'r')
semilogy(SNR_Vec, ber2)
grid
legend('Non-Equalized', 'Theoretical BER, No ISI','Equalized', 'location', 'southwest')
title(['BER Curves for ' num2str(m) '-ary QAM in Moderate ISI Channel']);
end
%% PART 2
clear all;

numIter = 1000;     %The number of iterations of the simulation
nSym = 1000;        %The number of symbols per packet
SNR = 12;           %12 dB SNR
M = 2;           

chan = [1 .2 .4];   %Moderate ISI


n = 31; %Codeword Length
k = 26; %Message Length
numBCH = 29;
training = 101; 

% Another set of values tested for BCH resulting in a higher BER
% n = 7; 
% k = 4; 
% numBCH = 132;
% training = 101;



% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, 1);

% eqobj = lineareq(6, rls(1)); % Attempt at Linear RLS

% Run the simulation numIter amount of times
for i = 1:numIter

    bits = randi([0,1], numBCH, k); %Generate random bits

    %BCH encoding
    bchmsg = gf(bits);
    code_a = bchenc(bchmsg,n,k);
    bits2 = reshape(code_a.x,1,numBCH*n);

%     msg = [randi([0 1], 1, training), bits2];
%     msg = reshape(msg,log2(M),nSym);
    msg = reshape([randi([0 1], 1, training), bits2],log2(M),nSym);
    
    tx = qammod(msg,M,'UnitAveragePower',true,'InputType','bit'); %BPSK modulate the signal
    
    txChan = filter(chan,1,tx); %Apply the channel
    txNoisy = awgn(txChan,SNR+10*log10(log2(M)), 'measured'); % Add AWGN

%   txNoisy2 = equalize(eqobj,txNoisy,tx(1:100)); %Attempt at Linerar Equalizer
                                                  % Note: Linear Equalizer gave
                                                  %a higher BER of about
                                                  %.5e-5
                                                 
    %Differential Feedback Equalizer
    stepsize = 0.01;
    dfeObj = dfe(5,3,lms(stepsize)); %Differential feedback object used for the Equalizer
    dfeObj.SigConst = qammod((0:M-1)',M,'UnitAveragePower',true)';
    dfeObj.ResetBeforeFiltering = 1;
    [tx_Eq,~,e] = equalize(dfeObj,txNoisy,tx(1:floor(training/(log2(M)))));


    rxMSG = qamdemod(tx_Eq,M,'UnitAveragePower',true,'OutputType','bit');%Demodulate

    %Decode
    channeled = reshape(rxMSG(training+1:end),numBCH,n);
    decoded = bchdec(gf(channeled),n,k);
    decbits = reshape(double(decoded.x),1,[]);

    %We're interested in the BER, which is the 2nd output of BITERR
    [~, berVec(i)] = biterr(bits(:), decbits(:));

end   %End numIter loop

%Compute the mean BER
ber = mean(berVec,1)
bitRate = k*numBCH/(nSym*log2(M))

%% Functions Used For Conversion
function [msg] = convertB2M(bits, M)
    % Convert message from bits into integer values.
    % M is a multiple of 2.
    length = log2(M); 
    msg = zeros(size(bits,1)/length, 1);  
    for i = 1:size(bits,1)/length
        %Using reccomended bi2de function
        msg(i) = bi2de(bits(1+(i-1)*length : 1+(i-1)*length + (length-1))');
    end
end
function [bits] = convertM2B(msg, M)
    % Convert message from integer values into bits.
    % M is a multiple of 2.
    length = log2(M); 
    bits = zeros(1, size(msg,1)*length);  
    for i = 1:size(msg,1)
        %Using reccomended bi2de function
        bits(1+(i-1)*length:1:1+(i-1)*length + (length-1)) = de2bi(msg(i), length);
    end
end
