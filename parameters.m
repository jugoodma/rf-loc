%%

% 1e9 GHz
% 1e6 MHz
% 1e3 KHz

% USRP N210 uses Gigabit Ethernet
% <https://kb.ettus.com/About_USRP_Bandwidths_and_Sampling_Rates#Host_Bandwidth>
% accordingly, we can stream 25 MS/s USRP <-> computer
% ~translates to 20 MHz of usable bandwidth
% this is baseband -- the RX chain is as follows:
%     RF frequency (detected by antenna)
%  -> complex baseband / IF signal
%  -> sampled by ADC
%  -> clocked into FPGA, down-conversion/decimation
%  -> streamed raw samples to computer

% we use 2.4 GHz antenna
% we use SBX-120 daughterboard
%  - has 400 MHz - 4.4 GHz freq coverage, 120MHz analog bandwidth
% we use USRP N210
%  - 100 MS/s ADC processing bandwidth
%  - 400 MS/s DAC processing bandwidth
% USRP N210 uses 1 Gigabit ethernet
%  - 25 MS/s @ 16-bit I/Q (full duplex)

% but essentially, the sampling rate is truly going to be
% number of samples collected / time spent collecting samples

% 20 MHz channel width
% 64 subcarriers
% -> 64/20e6 = 3.2e-6 seconds
% -> 1/3.2e-6 = 312.5e3 = 312.5 KHz frequency gap between subcarriers

% 2.4 GHz center frequency
% C = 3e8 m/s electromagnetic wave speed
% -> C/fc = 0.125 meters / wavelength
% -> 0.125 m/w * 100 cm/m /2 = 6.25 centimeters per half-wavelength
% == antenna distance

% ofdm data

bits_pilots = [1 0 1 0];
bits_preamble = [...
                            0 1 1 0 ...
    0 1 0 1 0 1 1 1 1 1 0 0 1 1 0 1 ...
    0 1 0 0 0 0 0 0 1 1 0 0 1 0 1 0 ...
    0 0 0 0 0 1 1 0 0 1 0 1 0 0 0 0 ...
];
cp = 16;
fc = 2.4e9;
fs = 20e6; % 10e6 given? USRP clock rate is 100MHz
guard_bins = [-32 -31 -30 -29 -28 -27 0 27 28 29 30 31];
num_bins = 64;
num_bins_data = 48;
num_bins_pilots = 4;
num_samples = 6288;
num_syms_data = 72;
num_syms_preamble = 8;
pilots = [-21 -7 7 21];
w = 312.5e3;

P = struct(...
    'bits_pilots',bits_pilots,...
    'bits_preamble',bits_preamble,...
    'cp',cp,...
    'fc',fc,...
    'fs',fs,...
    'guard_bins',guard_bins,...
    'num_bins',num_bins,...
    'num_bins_data',num_bins_data,...
    'num_bins_pilots',num_bins_pilots,...
    'num_samples',num_samples,...
    'num_syms_data',num_syms_data,...
    'num_syms_preamble',num_syms_preamble,...
    'pilots',pilots,...
    'w',w...
);


% general data

IntDeci_factor = 100;
SamplesPerFrame = 6500;
c = 3e8; % physconst('LightSpeed')
d = c/fc/2; % meters

