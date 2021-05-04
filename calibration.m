%%
% run this to test if the USRPs are calibrated/synchronized
%
% setup: have 1 USRP transmitter and 3 USRP receivers. we use USRP N210.
% the receiver radios should have antennas set up at half the
% center frequency (transmission) wave-length. For our set-up, we have
% $f_c = 2.4$GHz, thus with $c=3 \times 10^8$ we have antenna spacing 
% $\lambda/2 = c/f_c/2 = 0.0625$ meters.
%
% this code
% - sets up USRP devices
% - transmits a pure tone at $f_c = 2.4$GHz
% - plots received data
%

parameters;
addpath ofdm;

freq = fs/2/100;
% we want the transmission frequency to be 2.4 GHz
% actual transmitted = f_c + freq
% 2.4GHz = f_c + 100KHz
fc = 2.4e9 - freq;

% TODO ???

disp('loading.');
tic

IntDeci_factor = 100;
SamplesPerFrame = 6500;

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
disp('loaded.');

%%

sinewave = dsp.SineWave(1, freq);
sinewave.SampleRate = fs; % or is it ADCSampleRate?
sinewave.SamplesPerFrame = SamplesPerFrame; % div by 2*100
sinewave.OutputDataType = 'double';
sinewave.ComplexOutput = true;
signal = step(sinewave);
figure
plot(real(signal))
plot_fft(signal,fs);


%%
pause(3)
rxLog = dsp.SignalSink;
disp("Transmission Started");
start = clock;
% for the first 10 seconds, it'll transmit whatever it already was
% transmitting. then, it'll switch to what it should be (ish)
while etime(clock, start) <= 10
    tx(signal);
    % rx() returns M columns
    [data,len] = rx();
    rxLog(complex(data)) % casting -> avoid errs
end
disp("Transmission Stopped");


%%

data = rxLog.Buffer;
size(data)

%%
% small window

% offset = 1.6e7+3000;
% data = double(rxLog.Buffer);
% offset = 1e7;
offset = 1;
% window = 1*SamplesPerFrame;
window = 1000;
figure
hold on
for i=1:M
    plot(real(data(offset:offset+window,i)))
end
title('RX ---> TX')
legend(rxIP)
hold off

% for i=1:M
%     plot_fft(data(:,i),fs);
% end


%%
% whole thing

figure
hold on
for i=1:M
    plot(real(data(:,i)))
end
legend(rxIP)
hold off

