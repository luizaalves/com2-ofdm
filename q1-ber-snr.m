close all;
clear all;
clc;

%%
M = 4; % Q-PSK

nIteracoes = 10^3; 

nFFT        = 64;       % fft size
nDSC        = 48;       % number of data subcarriers
nSym        = 10^1;     % number of symbols per transmission
nPilot      = 4;        % number of pilot subcarriers
nPrefix     = 16;       % cyclic prefix

EbN0dB      = [-5:5:50]; % bit to noise ratio
EsN0dB      = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(nFFT/(nFFT+nPrefix)) % converting to symbol to noise ratio

BER = zeros(length(M),length(EbN0dB));
BER_M = zeros(length(M),length(EbN0dB));


    
nBitPerSym  = log2(M)*nDSC; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)


mod = comm.PSKModulator(M,'BitInput',true);
demod = comm.PSKDemodulator(M,'BitOutput',true);


%% Transmitter
message = randi([0 1],nSym*nBitPerSym,1);

% M-PSK Modulation
DataPSK = step(mod,message);

% Convers?o S?rie-Paralelo
modData = reshape(DataPSK,nDSC,nSym);

% OFDM subcarriers
SubCar = zeros(nDSC+nPilot+1,nSym); % Data + Pilot + DC component
% Add the pilot subcarriers:
SubCar(6,:) = 1;
SubCar(20,:) = 1;
SubCar(34,:) = 1;
SubCar(48,:) = -1;
% Add data to the subcarriers:
SubCar(1:5,:) = modData(1:5,:);
SubCar(7:19,:) = modData(6:18,:);
SubCar(21:26,:) = modData(19:24,:);
SubCar(27,:) = 0; % DC
SubCar(28:33,:) = modData(25:30,:);
SubCar(35:47,:) = modData(31:43,:);
SubCar(49:53,:) = modData(44:48,:);

OFDM_data = zeros(nFFT,nSym);
OFDM_data(1:6,:) = 0;        % zero PAD
OFDM_data(7:59,:) = SubCar;  % Data + Pilot + DC
OFDM_data(60:nFFT,:) = 0;    % zero PAD

% IFFT + Cyclic prefix
TX_OFDM = zeros(nFFT+nPrefix,nSym);

OFDM_ifft = ifft(OFDM_data);
TX_OFDM(1:nPrefix,:) = OFDM_ifft(nFFT-nPrefix+1:nFFT,:);
TX_OFDM(nPrefix+1:nFFT+nPrefix,:) = OFDM_ifft;

clear OFDM_ifft SubCar;
% Conversao Paralelo-Serie
TX = reshape(TX_OFDM,1,(nFFT+nPrefix)*nSym);
    
%% Channel    
for ni = 1:length(EbN0dB)
    for nInter = 1:nIteracoes

        % AWGN Channel
        Y = awgn(TX,EbN0dB(ni));

%% Receiver            
        % Conversao Serie-Paralelo
        R = reshape(Y,(nFFT+nPrefix),nSym);

        % OFDM Receiver
        % Remove the cyclic prefix
        RX_OFDM =  R(nPrefix+1:nFFT+nPrefix,:);

        % FFT
        RX = fft(RX_OFDM);
        RX_DataPiltotDC = RX(7:nFFT-5,:);    % Remove the zero PADs

        % Remove the Pilots and the DC component
        RX_data(1:5,:) = RX_DataPiltotDC(1:5,:);
        RX_data(6:18,:) = RX_DataPiltotDC(7:19,:);
        RX_data(19:24,:) = RX_DataPiltotDC(21:26,:);
        RX_data(25:30,:) = RX_DataPiltotDC(28:33,:);
        RX_data(31:43,:) = RX_DataPiltotDC(35:47,:);
        RX_data(44:48,:) = RX_DataPiltotDC(49:53,:);

        yMod = reshape(RX_data,nSym*nDSC,1);

        clear RX_DataPiltotDC;

        % M-PSK Demodulation
        yDemod = step(demod,yMod);

        indBER = find(message ~= yDemod);
        Ber_M = length(indBER);

        BER_M(1,ni) = BER_M(1,ni) + Ber_M;
    end

    BER(1,ni) = BER_M(1,ni)/(2*nSym*nBitPerSym*nIteracoes) 

end

close all; 
figure
semilogy(EbN0dB,BER(1,:),'bx-','LineWidth',2);
hold on;
grid on
xlabel('SNR, dB')
ylabel('Bit Error Rate')
title('BER for Q-PSK using OFDM')
legend('QPSK','Location','southwest');

