%NÃO FUNCIONOU
modOrder = 4;  % for 16-QAM
bitsPerSymbol = log2(modOrder)  % modOrder = 2^bitsPerSymbol
nFFT=64;
numCarr = 48;  % number of subcarriers
cycPrefLen = 16;  % cyclic prefix length
mpChan = [0.8; zeros(7,1); -0.5; zeros(7,1); 0.34];  % multipath channel
SNR = 15;
%SNR = [-5:5:50]; % bit to noise ratio
%SNR      = SNR + 10*log10(numCarr/nFFT) + 10*log10(nFFT/(nFFT+cycPrefLen)) % converting to symbol to noise ratio   % dB, signal-to-noise ratio of AWGN

nIteracoes = 10^3; 

numGBCarr = numCarr/16
gbLeft = 1:numGBCarr 
gbRight = (numCarr-numGBCarr+1):numCarr
dcIdx = (numCarr/2)+1
nullIdx = [gbLeft dcIdx gbRight]'
numDataCarr = numCarr - length(nullIdx)
numBits = numDataCarr*bitsPerSymbol;

srcBits = randi([0,1],numBits,1);
mod = comm.PSKModulator(modOrder,'BitInput',true);
qpskmod= step(mod,srcBits);
ofdmmod = comm.OFDMModulator('FFTLength',64, ...
    'PilotInputPort',true, ...
    'InsertDCNull',true, ...
    'CyclicPrefixLength',16, ...
    'NumSymbols', 10^1, ...
    'NumTransmitAntennas',2);
ofdmdemod = comm.OFDMDemodulator('FFTLength',64,'CyclicPrefixLength',16);

ofdmModOut = ofdmmod(qpskmod);

for ni = 1:length(SNR)
    for nInter = 1:nIteracoes
        mpChanOut = filter(mpChan,1,ofdmModOut);
        chanOut = awgn(mpChanOut,SNR(ni));

        ofdmDemodOut = ofdmdemod(chanOut,numCarr,cycPrefLen,cycPrefLen,nullIdx);
        
        mpChanFreq = fftshift(fft(mpChan,numCarr));
        mpChanFreq(nullIdx) = [];
        eqOut = ofdmDemodOut ./ mpChanFreq;
        scatterplot(eqOut)
        title('Frequency Domain Equalizer Output')

        demod = comm.PSKDemodulator(modOrder,'BitOutput',true);
        qpskdDemodOut=demod(eqOut);

        indBER = find(srcBits ~= qpskdDemodOut);
        Ber_M = length(indBER);

        Ber_M(1,ni) = Ber_M(1,ni) + Ber_M;

    end
  BER(1,ni) = BER_M(1,ni)/(2*numBits*nIteracoes) 
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
