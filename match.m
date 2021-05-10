%%
%
% goal: take in .dat files, which contain the actual time-domain signal,
% and return AoA/ToF estimates using spotfi
%
% assume samples are synchronized (recorded at same time, and radios are
% clock-synchronized)
%

addpath ofdm
addpath spotfi

parameters; % load variables
% loads P variable, which is a struct containing all OFDM variables
% loads other variables used during capture
% TODO

% size: samples X antennas
filename = "data/recording-45in-to-20in-150deg-to-90deg-moving-7sec.dat";
disp(strcat('Reading data file [',filename,'].'));
tic
samples = readmatrix(filename);
toc
M = 3; % width of samples, number of antennas

%
% re-calculate sent bits
%
rng(12345);
bits = zeros(1, num_syms_data*num_bins_data);
repeated = randi([0,1], num_bins_data, 1);
for i=1:num_syms_data
    bits(1+(i-1)*num_bins_data:i*num_bins_data) = repeated;
end
% this is the data packet per-symbol:
% repeated = [1 0 0 0 1 1 1 1 1 1 1 1 ...
%             0 0 0 1 1 1 1 1 1 1 0 0 ...
%             0 1 1 1 1 0 0 1 1 0 0 1 ...
%             1 0 1 1 1 1 1 0 0 0 0 0];

%%

parpool

%%
% do analysis
%

num_recorded_frames = length(samples)/SamplesPerFrame;
num_packets = 50; % pre-allocation length for arrays
packet_indices = zeros(M,num_packets); % indices in samples where 
%start_search_index = 2.5e6; % in case you have weird data at the start
start_search_index = 88e4;
%start_search_index = 950050;
index = ones(M,1)*start_search_index;

% spotfi expects H in this order
% SubCarrInd = [-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,-1,1,3,5,7,9,11,13,15,17,19,21,23,25,27,28]; % from IEEE 802.11n-2009 standard
SubCarrInd = [-26:-1 1:26]; % use all non-guard subcarriers
ind = fftshift(convert_bin_index_normal_to_fft(SubCarrInd,num_bins));
N = length(SubCarrInd);
CSI = zeros(N,M);
CSI_extended = zeros(N,M,num_syms_data);
bits_data = zeros(1,num_syms_data*num_bins_data);
BER  = ones(M,num_packets);
bidx = ones(M,1);
aoa = zeros(2,num_packets); % use 2 here b/c spotfi returns 2 values
tof = zeros(2,num_packets);
aoa_extended = zeros(2,num_syms_data,num_packets); % each packet has array of AoA/ToF
tof_extended = zeros(2,num_syms_data,num_packets);

window = ceil(num_samples*1.5); % search window of given sample
shift  = ceil(num_samples/2);
subcarrier_config = ones(1,num_bins);
subcarrier_config(convert_bin_index_normal_to_fft(guard_bins,num_bins)) = 0;
subcarrier_config(convert_bin_index_normal_to_fft(pilots,num_bins)) = 3;
done = 0;
j = 1;
examine = 1:M; % the antennas we need to search
alignmentThreshold = 100;

% spotfi unchanging variables
% https://dhalperi.github.io/linux-80211n-csitool/ -> 30 subcarriers for H
% https://wands.sg/research/wifi/AtherosCSI/ -> 52 subcarriers for H
% 20 MHz channel bandwidth / 64 total bins = 312.5 KHz freq gap
fgap = 312.5e3; % frequency gap in Hz between successive subcarriers in WiFi
lambda = c/fc; % wavelength = 2*d
T = 1; % number of transmitter antennas
% MUSIC algorithm requires estimating MUSIC spectrum in a grid. paramRange captures parameters for this grid
% For the following example, MUSIC spectrum is caluclated for 101 ToF (Time of flight) values spaced equally between -25 ns and 25 ns.
% MUSIC spectrum is calculated for for 101 AoA (Angle of Arrival) values between -90 and 90 degrees.
paramRange = struct;
paramRange.GridPts = [101 101 1]; % number of grid points in the format [number of grid points for ToF (Time of flight), number of grid points for angle of arrival (AoA), 1]
paramRange.delayRange = [-50 50]*1e-9; % lowest and highest values to consider for ToF grid. 
%paramRange.angleRange = 90*[-1 1]; % lowest and highest values to consider for AoA grid.
paramRange.angleRange = [0 360]; % as a circle
do_second_iter = 0;
% paramRange.seconditerGridPts = [1 51 21 21];
paramRange.K = floor(M/2)+1; % parameter related to smoothing.  
paramRange.L = floor(N/2);   % parameter related to smoothing.  
paramRange.T = 1;
paramRange.deltaRange = [0 0];
maxRapIters = Inf;
useNoise = 0;
paramRange.generateAtot = 1;
% 0 = load saved matrix
% 1 = don't save matrix
% 2 = save matrix

% MAIN LOOP
tic
while ~done
    % detect and decode packets
    % compute channel estimates
    for i=examine % antenna samples we need to check
        if i == 0; continue; end % skip zeros (set during alignment check)
        found = 0; % have not found a packet
        % detect packet
        while ~found && index(i) <= length(samples)-window % sliding window
            candidate = transpose(samples(index(i):index(i)+window,i));
            idx = packet_detection(candidate,P); % only need short window
            if idx == -1
                % no packet found
                index(i) = index(i) + shift;
                continue;
            end
            % packet found!
            found = 1; % set flag to leave detection loop
            packet_indices(i,j) = idx + index(i); % save packet index
            % estimate preamble channel (standard OFDM)
            candidate = transpose(samples(index(i):end,i));
            current_index = idx + 4*num_bins;
            cfo = estimate_cfo(...
                candidate(current_index:current_index+num_bins-1),...
                candidate(current_index+num_bins:current_index+num_bins+num_bins-1),P);
            candidate(current_index:end) = correct_cfo(candidate(current_index:end),cfo,P);
            current_index = current_index + 2*num_bins;
            H = estimate_channel(...
                candidate(current_index:current_index+num_bins-1),...
                candidate(current_index+num_bins:current_index+2*num_bins-1),P);
            CSI(:,i) = transpose(H(ind)); % save channel estimate H
            % compute channel throughout the packet, and decode bits (standard OFDM)
            % decoding does not handle indexing errors
            current_index = current_index + 2*num_bins+cp;
            for m = 1:num_syms_data % decode loop
                symt = candidate(current_index:current_index+num_bins-1);
                symf = (1/sqrt(num_bins))*fft(symt);
                [r_cfo,r_sfo] = estimate_residual_cfo_sfo(symf,H,P);
                H = correct_residual_cfo_sfo(H,r_cfo,r_sfo,P);
                % TODO -- we might know the sent data as well, can we
                % compute the channel that way? or is it possible that each
                % TX pedestrian sends some sort of unique identifier?
                CSI_extended(:,i,m) = transpose(H(ind)); % save channel estimate H from this data symbol
                symf = correct_channel(symf,H,P);
                decoded = symf(subcarrier_config==1) < 0;
                bits_data(1+(m-1)*num_bins_data:m*num_bins_data) = fftshift(decoded);
                current_index = current_index + num_bins+cp;
            end
            % compute Bit Error Rate
            cbits = sum(bits == bits_data);
            BER(i,bidx(i)) = 100*(1-(cbits/length(bits)));
            % keep searching
            index(i) = index(i) + current_index;
            % with skipping windows -- we can't skip num_samples ahead
            % because sometimes the packet isn't num_samples ahead
            % (ie, if the packet was broken in the middle, or if some
            % samples were skipped or something)
            % ---
        end % packet detection loop
    end % antenna loop
    % skip packets that yield non-zero BER
    % assert bidx() is all the same
    if sum(bidx == bidx(1)) < length(bidx)
        fprintf('b');
    end
    if sum(BER(:,bidx(1))) > 0 % threshold?
        fprintf('e');
        continue;
    end
    for i=1:M
        bidx(i) = bidx(i)+1;
    end
    % check alignment of packets
    %  we expect packets to be miss-aligned by only a few samples
    % for now, let's just require them to be within 100 samples
    % if they're not, we just keep searching for the
    % lagged-behind antennas (so just skip this packet for now)
    % IE, update the "examine" array with the arrays we need to catch up
    %
    % first calculate the maximum index, with corresponding antenna
    [sorted_indices,index_of_sort_index] = sort(index,'descend');
    sorted_indices = sorted_indices(1) - sorted_indices; % eg: [0, 10, 101]
    % any that are above threshold we must re-examine
    examine=zeros(1,M);
    for i=1:M
        if sorted_indices(i) > alignmentThreshold
            examine(i) = index_of_sort_index(i);
        end
    end
    % are there antennas that need to be caught up?
    if sum(examine) > 0
        for i=1:M
            if index(i) > length(samples)-window
                done = 1;
            end
        end
        % print out the problematic index set, and continue main loop
        fprintf(strcat("\n [",strjoin(string(index)),"]\n"));
        continue;
    end
    % else, detected packets are aligned, so we can compute spotfi estimates
    examine = 1:M; % reset examine array to all M antennas
    % spotfi
    % https://bitbucket.org/mkotaru/spotfimusicaoaestimation/src/master/
    CSI_reshaped = reshape(CSI,N*M,1); % reshape channel state information for spotfi
    % ToF sanitization code (Algorithm 1 in SpotFi paper)
    csi_plot = reshape(CSI_reshaped, N, M);
    [PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
    ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
    csi_plot = csi_plot.*ToMult;
    relChannel_noSlope = reshape(csi_plot, N, M, T);
    sample_csi_trace_sanitized = relChannel_noSlope(:);
    % MUSIC algorithm for estimating angle of arrival
    % aoaEstimateMatrix is (nComps x 5) matrix where nComps is the number of paths in the environment.
    % First column is ToF in ns and second column is AoA in degrees as defined in SpotFi paper
    aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized, M, N, c, fc,...
            T, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(2));
    % save estimates to our arrays
    tof(:,j) = aoaEstimateMatrix(:,1); % ToF in nanoseconds
    aoa(:,j) = aoaEstimateMatrix(:,2); % AoA in degrees
    
    % repeat spotfi for each channel estimate from all symbols
    % THIS TAKES A WHILE (but CAN BE PARALLELIZED ez)
    parfor m=1:num_syms_data
        CSI_reshaped = reshape(CSI_extended(:,:,m),N*M,1); % reshape channel state information for spotfi
        % ToF sanitization code (Algorithm 1 in SpotFi paper)
        csi_plot = reshape(CSI_reshaped, N, M);
        [PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
        ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
        csi_plot = csi_plot.*ToMult;
        relChannel_noSlope = reshape(csi_plot, N, M, T);
        sample_csi_trace_sanitized = relChannel_noSlope(:);
        % MUSIC algorithm for estimating angle of arrival
        % aoaEstimateMatrix is (nComps x 5) matrix where nComps is the number of paths in the environment.
        % First column is ToF in ns and second column is AoA in degrees as defined in SpotFi paper
        aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized, M, N, c, fc,...
                T, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(2));
        % save estimates to our arrays
        tof_extended(:,m,j) = aoaEstimateMatrix(:,1); % ToF in nanoseconds
        aoa_extended(:,m,j) = aoaEstimateMatrix(:,2); % AoA in degrees
    end

    % print out progress
    fprintf('.');
    if mod(j,50)==0
        fprintf('\n')
    end
    % leave if necessary (unsure if needed)
    for i=1:M
        if index(i) > length(samples)-window
            done = 1;
        end
    end
    j = j + 1;
end
toc
disp('Done.');


%%
% show found packets along with sampled data
%

%trim = start_search_index;
trim = 1;
figure('Renderer', 'painters', 'Position', [10 10 900 600])
for i=1:M
    subplot(M,1,i)
    hold on
    plot(trim:length(samples(:,i)),real(samples(trim:end,i)),'k')
    %plot(88e4:92e4,real(samples(88e4:92e4,i)))
    %plot(packet_indices(i,BER(i,:) == 0),zeros(1,length(BER(i,BER(i,:) == 0))),'r*')
    % end-1 in case the previous code errored on the last packet
    plot(packet_indices(i,1:end-1),real( samples(packet_indices(i,1:end-1),i) ),'r*')
    %title(strcat("RX Antenna ",num2str(i+1)," - ",filename))
    title(strcat("RX Antenna ",num2str(i)))
    xlabel('Samples')
    ylabel('Amplitude')
    xlim([1 length(samples(:,i))])
    hold off
end

%%
% view AoA estimates
%

% put angle into
%     ^
%     |
% 2   3   4
% 180deg  0deg

figure('Renderer', 'painters', 'Position', [10 10 900 600])
aoaidx = 1:length(aoa);
aoaidx = aoaidx(sum(BER)==0);
for i=aoaidx
    subplot(2,1,1)
    hold on
    plot(1+(i-1)*num_syms_data:i*num_syms_data,mod(transpose(aoa_extended(1,:,i))-90,180),'k.')
    plot(1+(i-1)*num_syms_data,mod(aoa(1,i)-90,180),'r*')
    %plot(1+(i-1)*num_syms_data:i*num_syms_data,mod(transpose(aoa_extended(1,:,i)),360),'k.')
    %plot(1+(i-1)*num_syms_data,mod(aoa(1,i),360),'r*')
    hold off
    
    subplot(2,1,2)
    hold on
    plot(1+(i-1)*num_syms_data:i*num_syms_data,transpose(tof_extended(1,:,i)),'k.')
    plot(1+(i-1)*num_syms_data,tof(1,i),'r*')
    hold off
end

subplot(2,1,1)
% legend({'whole packet','just preamble'})
title(strcat("Spotfi AoA Estimate"))
% axis([1 num_syms_data*length(aoa) -90 90])
axis([1 num_syms_data*length(aoa) 0 360])
ylabel('Degrees')
xlabel('Symbols')

subplot(2,1,2)
% legend({'whole packet','just preamble'})
title(strcat("Spotfi ToF Estimate"))
axis([1 num_syms_data*length(tof) -50 50])
ylabel('Nano-seconds')
xlabel('Symbols')

%%
% distance/angle from center antenna
% angle uses 90 deg == ray perpendicular to center antenna
%  antenna 3->2 == 0 degrees
%  antenna 3->4 == 180 degrees
%
%     ^
%     |
% 2   3   4
% 0deg    180deg

% ~34 in == 0.8636 m
% 0.8636 m / 3e8 m/s * 1e9 ns/s = 2.8787 nano-seconds

% recording-34in-85deg-stationary.dat
% AOA -1.8
% offset 2e6

% recording-34in-175deg-stationary.dat
% offset = 2.535e6-600;
% AOA maybe 28? -- highly inaccurate answers

% recording-34in-265deg-stationary.dat
% AOA -3.6

% i'd expect the 85 and 265 to be opposites of each other. i guess not
% though?

% recording-45in-67deg-stationary-laptop-obstacle.dat
% AOA -14.4000
% close enough.

% recording-45in-to-20in-150deg-to-90deg-moving-7sec.dat
% AOA range [-9 -14.4]
% ??

