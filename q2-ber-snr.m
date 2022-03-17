%não está funcionando
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
EsN0dB      = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(nFFT/(nFFT+nPrefix)); % converting to symbol to noise ratio

BER = zeros(length(M),length(EbN0dB));
BER_M = zeros(length(M),length(EbN0dB));


    
nBitPerSym  = log2(M)*nDSC; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)


mod = comm.QPSKModulator;
demod = comm.QPSKDemodulator;

ofdmMod = comm.OFDMModulator('FFTLength',64, ...
    'PilotInputPort',true, ...
    'InsertDCNull',true, ...
    'CyclicPrefixLength',16, ...
    'NumSymbols', 10^1, ...
    'NumTransmitAntennas',2);
ofdmDemod = comm.OFDMDemodulator(ofdmMod);
ofdmDemod.NumReceiveAntennas = 2;
ofdmModDim = info(ofdmMod);

numData = ofdmModDim.DataInputSize(1);   % Number of data subcarriers
numSym = ofdmModDim.DataInputSize(2);    % Number of OFDM symbols
numTxAnt = ofdmModDim.DataInputSize(3);  % Number of transmit antennas
nframes = 100;
data = randi([0 3],nframes*numData,numSym,numTxAnt);
modData = mod(data(:));
modData = reshape(modData,nframes*numData,numSym,numTxAnt);

for k = 1:nframes

     % Find row indices for kth OFDM frame
    indData = (k-1)*ofdmModDim.DataInputSize(1)+1:k*numData;

    % Generate random OFDM pilot symbols
    pilotData = complex(rand(ofdmModDim.PilotInputSize), ...
        rand(ofdmModDim.PilotInputSize));

    % Modulate QPSK symbols using OFDM
    dataOFDM = ofdmMod(modData(indData,:,:),pilotData);
    
    
    
    %% Channel    
    for ni = 1:length(EbN0dB)
        for nInter = 1:nIteracoes

            % AWGN Channel
            Y = awgn(dataOFDM,EbN0dB(ni));

%% Receiver            
        % Conversao Serie-Paralelo
        R = reshape(Y,(nFFT+nPrefix),nSym,2);

        % OFDM Receiver
        % Remove the cyclic prefix
        yMod =  ofdmDemod(dataOFDM);

        % Q-PSK Demodulation
        yDemod = demod(yMod(:));

        
        indBER = find(data ~= yDemod);
        Ber_M = length(indBER);

        BER_M(1,ni) = BER_M(1,ni) + Ber_M;
        end
    end

    BER(1,ni) = BER_M(1,ni)/(2*nSym*nBitPerSym*nIteracoes) 

end
