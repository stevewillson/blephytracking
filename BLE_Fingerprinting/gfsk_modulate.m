function outputSignal  = gfsk_modulate(bitPattern, freqSep, samplingRate)
% Generate a BLE signal.
% bitPattern: pattern to convert to the BLE signal
% freqSep: separation frequency of 0 and 1 bits.
% samplingRate: Sampling rate (Hz).

% calculate the number of samples per bit
nSamplesPerBit = samplingRate / 1e6;

% make a vector for the sample_times that is increasing from time 1 to the
% total number of samples
secondsPerSample = 1 / samplingRate;

% make a vector that is increasing from 1 to the total number of samples 
% with each element value that is equal to the time
sampleTimes = (1:(nSamplesPerBit*length(bitPattern))) * secondsPerSample;

% allocate space for a frequency shift keyed array
gammaFsk = zeros(1,length(sampleTimes));
for i=1:length(bitPattern)
    startIndex = ((i-1) * nSamplesPerBit) + 1;
    endIndex = i * nSamplesPerBit;
    % assign a portion of the matrix to either 0 or 1 based on the current
    % value of the bit pattern
    gammaFsk(startIndex:endIndex) = ((bitPattern(i)*2)-1);
end

% create a lowpass Finite Impulse Response (FIR) Gaussian pulse-shaping
% filter, gauss_filter is a vector of filter coefficients
gaussFilter = gaussdesign(0.3, 3, nSamplesPerBit);
% figure; plot(gaussFilter);


% filter the gamma_fsk signal (0's or 1's) by the gauss_filter coefficients
gammaGfsk = filter(gaussFilter, 1, gammaFsk);
% figure; plot(gammaGfsk);
% calculate the phase of the signal, use cumulative trapezoidal integration
gfskPhase = (freqSep / samplingRate) * pi * cumtrapz(gammaGfsk);
outputSignal = exp(1i * gfskPhase);
% figure; plot(outputSignal);

% transpose the matrix, convert from a row to a column
outputSignal = outputSignal.';
% figure; plot(outputSignal);

end