function [signalFingerprint, bits] = ble_fingerprint(signal, snr, samplingRate, preambleDetect, interpolationFactor, nPartition)
% Compute a fingerprint from a detected BLE signal
% including CFO, I/Q offset and I/Q imbalance with several different methods

% center frequency of receiver
receiverCenterFrequency = 2.402e9;

% BLE channel
bleChannelFrequency = 2.402e9;

% upsampling factor for the BLE signal
upsampledSamplingRate = samplingRate * interpolationFactor;
nSamplesPerBit = upsampledSamplingRate / 1e6;

% normalize the signal by the average amplitude
normalizedSignal = signal / mean(abs(signal));

% the MATLAB 'normalize' function centers the values at 0,0, do not use
% that funciton to normalize the signal

% interpolate the signal
interpolatedSignal = interp(normalizedSignal,interpolationFactor);

% Remove the channel offset and move the signal fft to center
if interpolationFactor ~= 1
    % shift the zero frequency component to the center of the spectrum
    signalFft = fftshift(fft(interpolatedSignal));
    lSignal = length(interpolatedSignal);

    signalFftCentered = zeros(lSignal,1);
    lCenter = floor(lSignal/2);
    lChannel = floor((bleChannelFrequency-receiverCenterFrequency)/upsampledSamplingRate*lSignal+lSignal/2);
    lBandwidth = floor((lSignal-1)/nSamplesPerBit);
    signalFftCentered(lCenter-lBandwidth:lCenter+lBandwidth) = signalFft(lChannel-lBandwidth:lChannel+lBandwidth);
    centeredSignal = ifft(ifftshift(signalFftCentered));
end

% Decode the BLE signal
[bleDecodeOutputSignal, signalFreq, bits] = ble_decoder(centeredSignal, upsampledSamplingRate, preambleDetect);

% plot the steps of the processed signal
% tiledlayout(2,2);
% nexttile();
% plot(signal);
% title("Original Signal");
% nexttile();
% plot(normalizedSignal);
% title("Normalized Signal");
% nexttile();
% plot(interpolatedSignal);
% title("Interpolated Signal");
% nexttile();
% plot(centeredSignal);
% title("Centered Signal");

% figure;
% plot(bleDecodeOutputSignal);
% title("BLE Decoded Output Signal");

% Estimate CFO for initialization using preamble averaging method
endIndex = nSamplesPerBit * 8;
preamble = signalFreq(1:endIndex-1);
estCfo = mean(preamble(1:end));

% if the estimated CFO is above a threshold, set it to 0
cfoThreshold = 100e3;
estCfo2 = estCfo;
if abs(estCfo2) > cfoThreshold
    estCfo2 = 0;
end

% Run the imperfection estimator code
% should there be another parameter ignored?
% the BLE_Imperfection_Estimator_NAGD has 13 return parameters
% TODO use a struct to pass information back from the
% BLE_Imperfection_Estimator_NAGD
[amp, epsilon, phi, I, Q, IQO, IQI, f0, phi_off, error, ~, ~, ~] = BLE_Imperfection_Estimator_NAGD(bleDecodeOutputSignal,bits,upsampledSamplingRate,estCfo2,0,0,0,0,1,snr,nPartition);

% TODO - what is the 'tt' variable
timePerSample = 1 / upsampledSamplingRate;
tt = 0:timePerSample:length(bleDecodeOutputSignal)*timePerSample-timePerSample;
estimatedSignal = bleDecodeOutputSignal .* exp(-1j*(2*pi*f0*tt'+phi_off/(360/(2*pi))));

% Check if the signal ellipse is within a range
try
    sig = estimatedSignal(randperm(length(estimatedSignal)));
    ell = fit_ellipse(-real(sig)/amp*5,3*imag(sig)/amp);
    flag = 1;
catch
    warning('Ill ellipse');
    flag = 0;
end


angsig = angle(estimatedSignal);
spl = 8;
quar = zeros(1,spl);
for sp = 1:spl/2
    e1 = (((angsig>((sp-1)*2*pi/spl))+(angsig<(sp*2*pi/spl)))==2);
    quar(1,sp) = mean(estimatedSignal(e1));
    e1 = (((angsig>(-pi+(sp-1)*2*pi/spl))+(angsig<(-pi+sp*2*pi/spl)))==2);
    quar(1,sp+spl/2) = mean(estimatedSignal(e1));
end

signalFingerprint = struct();

signalFingerprint.error = error;
signalFingerprint.amp = amp;
signalFingerprint.f0 = f0;
signalFingerprint.estCfo = estCfo;
signalFingerprint.IQO = IQO;
signalFingerprint.I = I;
signalFingerprint.Q = Q;
signalFingerprint.iqMagnitude = sqrt(I^2+Q^2);
signalFingerprint.IQI = IQI;
signalFingerprint.epsilon = epsilon;
signalFingerprint.phi = phi;
signalFingerprint.filename = '';


% TODO add a fingerprint structure
% fingerprint_vec = [error,...
%                     amp,...
%                     f0,...
%                     est_cfo,...
%                     IQO,...
%                     I,...
%                     Q,...
%                     sqrt(I^2+Q^2),...
%                     IQI,...
%                     epsilon,...
%                     phi,...
%                     ell.X0/ell.a,...
%                     ell.Y0/ell.b,...
%                     ell.X0_in/ell.a,...
%                     ell.Y0_in/ell.b,...
%                     sqrt((ell.X0/ell.a)^2+(ell.Y0/ell.b)^2),...
%                     sqrt((ell.X0_in/ell.a)^2+(ell.Y0_in/ell.b)^2),...
%                     ell.a*3/ell.b/5,...
%                     ell.phi,...
%                     real(mean(quar)),...
%                     imag(mean(quar)),...
%                     abs(mean(quar)),...
%                     mean(real(signal)),...
%                     mean(imag(signal)),...
%                     abs(mean(real(signal))+1i*mean(imag(signal)))];
%
% % check if the length of the fingerprint vector is not 25, if not, set the
% % flag to 0
% if length(fingerprint_vec) ~= 25
%     flag=0;
% end
%
%
% % if the fingerprint vector was produced and the error is less than a
% % threshold, then set the fingerprint for return
% % error_threshold = 0.45;
% error_threshold = 0.5;
% if error(end) < error_threshold && flag == 1
%     fingerprint = [fingerprint;fingerprint_vec];
% end
