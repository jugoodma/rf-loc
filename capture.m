%%
% this file captures data from the radio setup
%
% setup:
% - 1 USRP N210 transmitter
% - 3 USRP N210 receivers
% - 1 linear antenna array properly spaced
% - 1 OctoClock-G CDA-2990
% -- this synchronizes the RX radio samples
% - 1 GiB network switch
%
% if you need to collect multiple sets of data, re-run this entire file
% top-bottom.
% you can also just release the tx/rx, then only re-run the
% recording/saving part.
%

% total transmission time
% NOTE -- the first 10 seconds of transmission are essentially null
% this is due to the drivers or something
% essentially, I think the computer has to establish a connection to the
% radios, and that takes a while.
% then, the radios have to start up their local oscillators or raise up
% their analog gain, or something
% then _eventually_ the signal is fine.
% we recommend using >= 30 seconds here.
seconds = input('How many seconds would you like to transmit? [30]: ');
if isempty(seconds)
    seconds = 30;
end

% this is the amount of pause time before transmission actually starts
% the script will tell you when it's pausing, and when it'll start
% transmitting
wait_time = input('How many seconds would you like to pause? [0]: ');
if isempty(wait_time) || wait_time < 0
    wait_time = 0;
end

filename = input('File name for recording (without extension) [rec.dat]: ','s');
if isempty(filename)
    filename = "rec";
end
filename = filename+".dat";


%
% get radios
%

disp('Instantiating radio objects.');
tic

addpath ofdm
addpath spotfi

parameters;

txIP = '192.168.10.5';
rxIP = [...
    "192.168.10.2",...
    "192.168.10.3",...
    "192.168.10.4"...
];
M = length(rxIP);
tx = comm.SDRuTransmitter(...
    'Platform','N200/N210/USRP2',...
    'IPAddress',txIP,...
    'CenterFrequency',fc,...
    'Gain',31.5,...
    'InterpolationFactor',IntDeci_factor...
);
rx = comm.SDRuReceiver(...
    'Platform','N200/N210/USRP2',...
    'IPAddress',convertStringsToChars(strjoin(rxIP,',')),...
    'ChannelMapping',1:M,...
    'CenterFrequency',fc,...
    'Gain',31.5,...
    'PPSSource','External',...    
    'ClockSource','External',...
    'DecimationFactor',IntDeci_factor,...
    'SamplesPerFrame',SamplesPerFrame...
);

toc
disp('Radios ready to go.');


%
% generate OFDM signal
% add a VERY small amount of space after the packet
% to fill up the entire frame
%
% TODO -- probably should be sending a continuous packet ping
% we have this here essentially for debugging
%

disp('Generating OFDM signal.');
rng(12345);
bits = zeros(1, num_syms_data*num_bins_data);
repeated = randi([0,1], num_bins_data, 1);
disp(strcat(...
    'repeating bits [',...
    transpose(num2str(repeated)) ,...
    '] for [',...
    num2str(num_syms_data),...
    '] symbols.'...
));
for i=1:num_syms_data
    bits(1+(i-1)*num_bins_data:i*num_bins_data) = repeated;
end
signal = tx_ofdm_chain(bits,P);
signal = signal / max(signal);

figure
plot(real(signal),'k');
axis([1 length(signal) -1 1])
title('Time domain OFDM packet');
xlabel('Samples')
ylabel('Amplitude')
plot_fft(signal, fs);
title('Frequency domain OFDM packet');

z = zeros(1, SamplesPerFrame);
z(1) = 1i;
signal = [signal zeros(1,SamplesPerFrame-length(signal))];
signal = transpose(signal);


%%
% actually transmit and receive the signal
%

disp('Ready to transmit/receive. Pausing before transmission.');
pause(wait_time)
rxLog = dsp.SignalSink;
disp('Transmission starting.');
start = clock;
% for the first 10 seconds, it'll transmit whatever it already was
% transmitting. then, it'll switch to what it should be (ish)
while etime(clock, start) <= seconds
    tx(signal); % pauses for a few seconds on first iteration
    [data,len] = rx(); % rx() returns M columns
    rxLog(complex(data)) % ignore len==0, save it all
    %fprintf('.'); % uncomment if you need to see when tx/rx actually starts
end
fprintf('\n');
disp('Transmission stopped.');
% using 3 RX, we have
% 13279500 samples / 15 seconds = 885300 Hz
% which is about 2.66 MHz, but it should be 25-ish MHz...
% weird. i still don't know why this happens.
% we use decimation factor == 100.
% I _think_ that
% NUM_SAMPLES * DECIMATION_FACTOR / SAMPLING_FREQUENCY / NUM_RX
% == total recording time
% i don't think this is totally right, but we're getting closer and closer

%
% save recording to disk
% release radios
%

% todo, is it faster to write the rxLog.Buffer, or to pull it out and
% _then_ write it...
disp('Writing data to file.');
tic
writematrix(rxLog.Buffer, filename);
disp('writing file took:');
toc

% to ensure that we don't get "previous data" during the next recording
disp('Releasing TX/RX.');
release(rx);
release(tx);

disp('Finished.');

