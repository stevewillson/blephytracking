function [amp, e, phi, I, Q, IQO, IQI, f0, phiOff, error, estSignal, fr, eO] = BLE_Imperfection_Estimator_NAGD(inputSignal, bits, samplingRate, initF0, initE, initPhi, initI, initQ, initAmp, signalToNoiseRatio, nPartitions)
% Use a Nesterov Accelerated Gradient Descent to estimate hardware
% imperfection fingerprints of a Bluetooth signal
% Source: "Evaluating Physical-Layer BLE Location Tracking Attacks on Mobile Devices", IEEE SP 22

% update snr to set the error threshold
updatedSignalToNoiseRatio = 10^(signalToNoiseRatio/20);
errorThreshold = max(0.4,1/(updatedSignalToNoiseRatio+1));

% x is the input signal
nSamplesPerBit = samplingRate / 1e6;

timePerSample = 1 / samplingRate;
lInputSignal = length(inputSignal);

inputSignal2 = inputSignal(1:end-(2*nSamplesPerBit));
t2 = (1:lInputSignal) * timePerSample;

% add 2 bits (1 0) to the beginning of the bits vector
bits = [1;0;bits];

% Create the clean signal
estPerfectSignal2 = gfsk_modulate(bits, 500e3, samplingRate).';
if nPartitions ~= 1
    estPerfectSignal2 = estPerfectSignal2((3.5 * nSamplesPerBit) + 1 : end - (0.5 * nSamplesPerBit));
else
    estPerfectSignal2 = estPerfectSignal2((3.0 * nSamplesPerBit) + 1 : end - (1 * nSamplesPerBit));
end

fr = [];
I2 = 0;
Q2 = 0;
IQO2 = 0;
IQI2 = 0;
e2 = 0;
phi2 = 0;
phiOff2 = 0;
f02 = 0;
amp2 = 0;
eO = [];
%freqsep2 = 0;
error2 = 0;
err = 0;

% Partition the samples to speed up the code and run over multiple
% random partitions to improve estimation accuracy and robustness
n = min(10, nPartitions);
l = floor(length(inputSignal2) / nPartitions);

for iPartition = 1:n
    % generate a random permutation of integers from 1 to length of
    % estPerfectSignal2 without repeating integers
    r = randperm(length(estPerfectSignal2));

    % only get the first l values (divided by the nPartitions)
    r = r(1:l);
    % assign estPerfectSignal to the first r samples of estPerfectSignal2
    estPerfectSignal = estPerfectSignal2(r);
    t = t2(r);
    inputSignal = inputSignal2(r);

    % initialization
    eNew = initE;
    phiNew = initPhi * pi / 180;
    iNew = initI;
    qNew = initQ;

    if iPartition == 1
        f0New = initF0;
    else
        f0New = f0;
        initF0 = f0;
    end
    w0New = 2*pi*f0New;
    phiOffsetNew = 36/360*2*pi;
    ampNew = initAmp;

    % Parameter update
    e = eNew;
    phi = phiNew;
    I = iNew;
    Q = qNew;
    w0 = w0New;
    phiOff = phiOffsetNew;
    amp = ampNew;
    %freqsep = freqsep_new;
    estSignal = ((amp-e)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+...
        1i*(amp+e)*(imag(estPerfectSignal)*cos(phi)+real(estPerfectSignal)*sin(phi))+...
        I+1i*Q) .* exp(1i*(w0*t+phiOff));

    eT = 0;
    phiT = 0;
    IT = 0;
    QT = 0;
    w0T = 0;
    phiOffT = 0;
    ampT = 0;

    errorDiff = 1;

    count = 0;
    round = 1;
    error = 1;

    while errorDiff > 1e-7 && count < 10e3 %2e3
        count = count + 1;

        % Set the learning rate and momentum
        if errorDiff > 1e-7
            learningRate = 1e-3;
            momentum = 0.9;

        elseif errorDiff < 1e-7
            learningRate = 1e-4;
            momentum = 0.9;

        else
            learningRate = 1e-3;
            momentum = 0.9;

        end

        % Momentum step
        %e = e_new - mom*e_t;
        phi = phiNew - (momentum * phiT);
        I = iNew - (momentum * IT);
        Q = qNew - (momentum * QT);
        w0 = w0New - (momentum * w0T);
        phiOff = phiOffsetNew - (momentum * phiOffT);
        %amp = amp_new - mom*amp_t;
        %freqsep = freqsep_new - mom*freqsep_t;

        % Re-initialization
        if count > round*2e2 && error(end) > errorThreshold
            if floor(round/2) * 2 == round - 1
                round = round + 1;
                f0 = initF0 - (floor(round/2) * 1.5e3);
                w0New = 2 * pi * f0;
            else
                round = round + 1;
                f0 = initF0 + (round/2) * 1.5e3;
                w0New = 2 * pi * f0;

            end
        end

        % Compute updated created imperfect signal
        ImagPart = ((amp-e)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+I).*sin(w0*t+phiOff) +...
            ((amp+e)*(imag(estPerfectSignal)*cos(phi)+...
            real(estPerfectSignal)*sin(phi))+Q).*cos(w0*t+phiOff);

        RealPart = ((amp-e)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+I).*cos(w0*t+phiOff) -...
            ((amp+e)*(imag(estPerfectSignal)*cos(phi)+...
            real(estPerfectSignal)*sin(phi))+Q).*sin(w0*t+phiOff);

        % Gradient descent calculation and update

        % freqsep_d = ((amp-e)*(*cos(phi)+imag(est_signal_perfect)*sin(phi))+I).*sin(w0*t+phi_off) +...
        %         ((amp+e)*(imag(est_signal_perfect)*cos(phi)+...
        %         real(est_signal_perfect)*sin(phi))+Q).*cos(w0*t+phi_off);

        eD =  -mean((imag(inputSignal.')-ImagPart).*((-(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))).*sin(w0*t+phiOff)+...
            (imag(estPerfectSignal)*cos(phi) +...
            real(estPerfectSignal)*sin(phi)).*cos(w0*t+phiOff)));

        phiD = -mean((imag(inputSignal.')-ImagPart).*((amp-e)*(-real(estPerfectSignal)*sin(phi)+imag(estPerfectSignal)*cos(phi)).*sin(w0*t+phiOff)+...
            ((amp+e)*(real(estPerfectSignal)*cos(phi) -...
            imag(estPerfectSignal)*sin(phi))).*cos(w0*t+phiOff)));

        ID = -mean((imag(inputSignal.')-ImagPart).*(sin(w0*t+phiOff)));

        QD = -mean((imag(inputSignal.')-ImagPart).*(cos(w0*t+phiOff)));


        w0D = -mean((imag(inputSignal.')-ImagPart).*(t.*((amp-e)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+I).*cos(w0*t+phiOff) - ...
            t.*((amp+e)*(imag(estPerfectSignal)*cos(phi)+...
            real(estPerfectSignal)*sin(phi))+Q).*sin(w0*t+phiOff)));

        phiOffD = -mean((imag(inputSignal.')-ImagPart).*(((amp-e)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+I).*cos(w0*t+phiOff) - ...
            ((amp+e)*(imag(estPerfectSignal)*cos(phi)+...
            real(estPerfectSignal)*sin(phi))+Q).*sin(w0*t+phiOff)));

        ampD = -mean((imag(inputSignal.')-ImagPart).* (((1)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+I).*sin(w0*t+phiOff) +...
            ((1)*(imag(estPerfectSignal)*cos(phi)+...
            real(estPerfectSignal)*sin(phi))+Q).*cos(w0*t+phiOff)));


        eD = eD - mean((real(inputSignal.')-RealPart).*((-(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))).*cos(w0*t+phiOff)-...
            (imag(estPerfectSignal)*cos(phi) +...
            real(estPerfectSignal)*sin(phi)).*sin(w0*t+phiOff)));

        phiD = phiD - mean((real(inputSignal.')-RealPart).*((amp-e)*(-real(estPerfectSignal)*sin(phi)+imag(estPerfectSignal)*cos(phi)).*cos(w0*t+phiOff)-...
            ((amp+e)*(real(estPerfectSignal)*cos(phi) -...
            imag(estPerfectSignal)*sin(phi))).*sin(w0*t+phiOff)));

        ID = ID - mean((real(inputSignal.')-RealPart).*(cos(w0*t+phiOff)));

        QD = QD + mean((real(inputSignal.')-RealPart).*(sin(w0*t+phiOff)));


        w0D = w0D - mean((real(inputSignal.')-RealPart).*(-t.*((amp-e)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+I).*sin(w0*t+phiOff) - ...
            t.*((amp+e)*(imag(estPerfectSignal)*cos(phi)+...
            real(estPerfectSignal)*sin(phi))+Q).*cos(w0*t+phiOff)));

        phiOffD = phiOffD - mean((real(inputSignal.')-RealPart).*(-((amp-e)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+I).*sin(w0*t+phiOff) - ...
            ((amp+e)*(imag(estPerfectSignal)*cos(phi)+...
            real(estPerfectSignal)*sin(phi))+Q).*cos(w0*t+phiOff)));

        ampD = ampD - mean((real(inputSignal.')-RealPart).* (((1)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+I).*cos(w0*t+phiOff) -...
            ((1)*(imag(estPerfectSignal)*cos(phi)+...
            real(estPerfectSignal)*sin(phi))+Q).*sin(w0*t+phiOff)));

        if errorDiff >1e-7
            eT = (momentum * eT) + (learningRate * eD);
            phiT = (momentum * phiT) + (learningRate * phiD);
            IT = (momentum * IT) + (learningRate/1 * ID);
            QT = (momentum * QT) + (learningRate/1 * QD);
            w0T = (momentum * w0T) + (1e8 * learningRate * w0D);
            phiOffT = (momentum * phiOffT) + (10 * learningRate * phiOffD);
            ampT = (momentum * ampT) + (learningRate * ampD);
            %freqsep_t = (momentum * freqsep_t) + (learningRate * freqsep_d);
        else
            eT = (momentum * eT) + (learningRate * eD);
            phiT = (momentum * phiT) + (learningRate * phiD);
            IT = (momentum * IT) + (learningRate/1 * ID);
            QT = (momentum * QT) + (learningRate/1 * QD);
            %w0_t = mom*w0_t + 1e8*lr*w0_d;
            phiOffT = (momentum * phiOffT) + (10 * learningRate * phiOffD);
            ampT = (momentum * ampT) + (learningRate * ampD);
        end

        %amp_new = amp + lr * mean((imag(x.')-Imag_part).* Imag_part / amp);

        eNew = eNew - eT;
        phiNew = phiNew - phiT;
        iNew = iNew - IT;
        qNew = qNew - QT;
        w0New = w0New - w0T;
        phiOffsetNew = phiOffsetNew - phiOffT;
        ampNew = ampNew - ampT;
        %freqsep_new = freqsep_new-freqsep_t;

        % update parameters
        e = eNew;
        phi = phiNew;
        I = iNew;
        Q = qNew;
        w0 = w0New;
        phiOff = phiOffsetNew;
        amp = ampNew;
        %freqsep = freqsep_new;

        estSignal = ((amp-e)*(real(estPerfectSignal)*cos(phi)-imag(estPerfectSignal)*sin(phi))+...
            1i*(amp+e)*(imag(estPerfectSignal)*cos(phi)+real(estPerfectSignal)*sin(phi))+...
            I+1i*Q) .* exp(1i*(w0*t+phiOff));

        % compute the error
        %error = [error,mean(abs(est_signal.' - x)./abs(x))/2];
        error = [error,mean(norm(estSignal.'-inputSignal).^2)/mean(norm(inputSignal).^2)/2];

        if length(error) > 1 && error(end) < errorThreshold
            errorDiff = abs(error(end)-error(end-1));
        end
    end
    err = max([err,error(end)]);

    signal = inputSignal.*exp(-1j*(w0*t'+phiOff));

    try
        ell = fit_ellipse(real(signal),3*imag(signal));
        IQO = sqrt((ell.X0/ell.a)^2+(ell.Y0/ell.b)^2);
        IQI = ell.a/ell.b*3;
        flag = 1;
    catch
        warning('Ill ellipse');
        flag = 0;
    end

    if flag == 1
        f0 = w0/(2*pi);
        phiOff = phiOff*(360/(2*pi));
        phi = phi *(360/(2*pi));
        fr = [fr;f0];
        I2 = I2+I;
        Q2= Q2+Q;
        IQO2 = IQO2+IQO;
        IQI2 = IQI2+IQI;
        %e2 = e2+e;
        e2 = e2+e;
        phiOff2 = phiOff2+phiOff;
        phi2 = phi2+phi;
        f02 = f02+f0;
        amp2 = amp2+amp;
        error2 = error2+err;
    else
        disp('Ellipse was not found.');
    end
end
I = I2/n;
I = -I/amp;
Q = Q2/n;
Q = Q/amp;
IQO = IQO2/n;
IQI = IQI2/n;
phi = phi2/n;
e = e2/n;
phiOff = phiOff2/n;
f0 = f02/n;
amp = amp2/n;
error = err;
