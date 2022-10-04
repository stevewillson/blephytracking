function [bleSignal, signalFreq, bits] = ble_decoder(signal, samplingRate, preambleDetect)
% Simple Bluetooth demodulator.
% More robust decoders can be used instead.
% The code assumes the input signal only consists of 1 Bluetooth signal.
% samplingRate: sampling frequency (Hz)
% preambleDetect: If set to 1, the decoder also decodes the preamble. If
% set to zero, it assumes the Bluetooth signal starts from the first sample

% freqSep - separation of 0 and 1 bits for gfsk signal
freqSep = 500e3;

% divide the samplingRate by 1M, calculate how many samples are
% used to specify a bit
nSamplesPerBit = samplingRate / 1e6;

% create preamble
% create an array of 11 elements that specify a BLE preamble
preambleBits = [0,1,0,1,0,1,0,1,0,1,0];

% create a GFSK modulated signal
preambleSignal = gfsk_modulate(preambleBits, freqSep, samplingRate);

startIndex = (nSamplesPerBit) * 2.5;
endIndex = length(preambleSignal) - (nSamplesPerBit) * 0.5;
% why only collect these samples when decoding the packet?
preambleSignal = preambleSignal(startIndex:endIndex);

% calculate the instantaneous frequency of the preamble
signalAngle = unwrap(angle(preambleSignal));
slope = signalAngle(2:length(signalAngle)) - signalAngle(1:length(signalAngle)-1);
preambleFreq = samplingRate * slope/(2*pi);
% figure; plot(preamble_freq);

% obtain the instantaneous frequency of the signal
signalAngle = unwrap(angle(signal));
slope = signalAngle(3:length(signalAngle)) - signalAngle(2:length(signalAngle)-1);
signalFreq = slope / (2*pi) * samplingRate;
signalFreq = [signalFreq;0];
% figure; plot(signal_freq);

% Find the preamble. 
% Assume the beginning of the packet has been detected almost
% accurately.
if preambleDetect == 0
    startIndex = 1;
else
    % find the length of the starting signal
    l = length(signalFreq);
    z = xcorr(signalFreq, preambleFreq);
    z = z(l+1:end);
    if length(z) > 20e-6 * samplingRate
        [~,startIndex] = max(abs(z(floor(2e-6 * samplingRate):floor(20e-6 * samplingRate))));
        %start_ind = start_ind + sampling_rate / 1e6 / 2;
    else
        startIndex = 1;
    end
end

    signal = signal(startIndex:end);
    signalFreq = signalFreq(startIndex:end);
    
    signalFreq = signalFreq(1:floor(length(signalFreq)/(nSamplesPerBit))*(nSamplesPerBit));
    bleSignal = signal(1:floor(length(signalFreq)/(nSamplesPerBit))*(nSamplesPerBit));

    nBits = length(signalFreq) / nSamplesPerBit;
    
    % create a matrix with a number of rows equal to the number of
    % bits and number of columns equal to the samples per bit
    bitsFreq = reshape(signalFreq, nSamplesPerBit, nBits)';

    % decode the value of the bits in the signal
    % for all samples in the bitsFreq matrix, take the average of the
    % middle two elements, if they are >0 then this represents a 1 bit
    % if they are less than 0, this represents a 0 bit
    bits = (mean(bitsFreq(:,nSamplesPerBit/2:nSamplesPerBit/2+1),2)>0);
    %bits = (mean(bits_freq(:,40:60),2)>0);
    %bits = (mode(bits_freq(:,1:100)>0,2)>0);
    
end