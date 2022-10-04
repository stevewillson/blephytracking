% This example reads 4 samples for CONNECT_REQ between a phone and a lock
% it estimates their hardware imperfection values using the alorithm
% described in the Evaluating Physical-Layer BLE Location Tracking Attacks
% on Mobile Devices paper.


% close all the opened figures
close all;

% clear the workspace
clear;

% parameter setup
% sample rate of 4MHz
sampleRate = 4e6;

signalToNoiseRatio = 40;
preambleDetect = 1;

% is the goal to shoot for a 100Ghz signal?
% originally, this was 32 with a 3.125e6 sampling frequency, I am sampling
% at 4MHz so change this to 25
% interp_fac = 32;
interpolationFactor = 25;

% what is the purpose of setting the number of partitions?
nPartitions = 250;
fingerprintSize = 25;

tic

inputFileList = ls("/home/user/test/saved_samples/phone*");
inputFileList = split(inputFileList);
nFiles = length(inputFileList) - 1;


% set up the struct
signalFingerprint = struct();

signalFingerprint.error = 0;
signalFingerprint.amp = 0;
signalFingerprint.f0 = 0;
signalFingerprint.estCfo = 0;
signalFingerprint.IQO = 0;
signalFingerprint.I = 0;
signalFingerprint.Q = 0;
signalFingerprint.iqMagnitude = 0;
signalFingerprint.IQI = 0;
signalFingerprint.epsilon = 0;
signalFingerprint.phi = 0;
signalFingerprint.filename = '';

% create an array of structs for the files that will be read
fingerprintAll = repmat(signalFingerprint,nFiles,1);

% fingerprint_all = zeros(nFiles,f1ingerprint_size);
% store the structs in this array
% fingerprint_all = zeros(nFiles,1);

% fingerprint_all = ();

inputFile = "";
for iFile = 1:nFiles
    % Read the file including the signal
    inputFile = inputFileList{iFile};
    
    % file = "/home/user/test/saved_samples/phone2_test-1662762488_269973"

    % open the sample file, set the file to fid
    fid = fopen(inputFile, 'r');

    % create a signal matrix, specify the precision as 'int' (8bits of data)
    % reads in the file starting at the beginning
    [signal, ~] = fread(fid, 'int8');

    % close the file handle
    fclose(fid);

    % moves samples 1,3,5,7... to the 2nd column
    signal = reshape(signal, 2, []).';

    % make the 2nd column imaginary numbers and add it to the first column
    % make the matrix into a column matrix
    complexSignal = signal(:,1) + 1i * signal(:,2);
    % : is shorthand for 1:end
    % extract the row 1 - 10 with all columns
    % signal(1:10,:);
    
    % this was end-12 before, just making it 'end'
    % TODO - ???

    % the last part of the signal may contain huge values that distort the
    % average value of the signal

    % the signal is only 8 bit resolution, does it need to be divided by a
    % constant factor?
    %     signal = signal(1:end);

    % plot the signal in the complex plane
    %     plot(signal, "red")

    % the plots seem to show that the first ~20 samples or so are not valid
    % parts of the signal (according to the plot of the signal in the
    % complex plane
    complexSignal = complexSignal(20:end);
    %     plot(signal)

    % signal normalization is done in the BLE_Fingerprint.m module

    % Physical layer fingerprinting
    [fingerprintAll(iFile), bits] = ble_fingerprint(complexSignal, signalToNoiseRatio, sampleRate, preambleDetect, interpolationFactor, nPartitions);

    % extract the filename from the full file path
    [~, filename, ~] = fileparts(inputFile);
    fingerprintAll(iFile).filename = filename;
end
toc

% plot the fingerprint values
nFingerprints = length(fingerprintAll);
nSamplesPerPhone = 4;
phone1_cfo = zeros(nSamplesPerPhone, 1);
phone1_iqOffset = zeros(nSamplesPerPhone, 1);
phone2_cfo = zeros(nSamplesPerPhone, 1);
phone2_iqOffset = zeros(nSamplesPerPhone, 1);

iPhone1 = 1;
iPhone2 = 1;

for iFingerprint = 1:nFingerprints
    estCfo = fingerprintAll(iFingerprint).estCfo;
    iqOffset = fingerprintAll(iFingerprint).IQO;
    filename = fingerprintAll(iFingerprint).filename;

    if (contains(filename, "phone2"))
    % the iPhone SE
        phone2_cfo(iPhone2) = estCfo;
        phone2_iqOffset(iPhone2) = iqOffset;
        iPhone2 = iPhone2 + 1;
    else
        % the Galaxy S8 Phone
        phone1_cfo(iPhone1) = estCfo;
        phone1_iqOffset(iPhone1) = iqOffset;
        iPhone1 = iPhone1 + 1;
    end
end

sPhone1 = scatter(phone1_cfo, phone1_iqOffset, 'diamond', 'filled'); hold on;
sPhone1.DisplayName = "phone1";
sPhone2 = scatter(phone2_cfo, phone2_iqOffset, 'o', 'filled'); hold on;
sPhone2.DisplayName = "phone2";
legend();
title("Plot of Estimated CFO vs IQ Offset")
xlabel("Estimated CFO");
ylabel("IQ Offset")
hold off;