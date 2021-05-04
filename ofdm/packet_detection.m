% function packet_start_index = packet_detection(rx_signal,P)
%     num_bins = P.num_bins;
%     cp = P.cp;
%     %num_syms_preamble = P.num_syms_preamble;
%     num_syms_preamble = min(5,P.num_syms_preamble);
%     guard_bins = P.guard_bins;
%     pilots = P.pilots;
%     bits_preamble = P.bits_preamble;
%     num_samples = P.num_samples;
%     % Dc removal filter
%     alpha = 0.975;
%     rx_signal = filter([1 -1], [1 -alpha], rx_signal);
%     threshold=10; % PeakThreshold = 0.6
%     % generate preamble
%     signal = zeros(1,num_syms_preamble*num_bins);
%     %signal = zeros(1,num_bins);
%     symbol_freq = zeros(1,num_bins);
%     subcarrier_config = ones(1,num_bins);
%     subcarrier_config(convert_bin_index_normal_to_fft(guard_bins,num_bins)) = 0;
%     subcarrier_config(convert_bin_index_normal_to_fft(pilots,num_bins)) = 3; 
%     for m = 1:num_syms_preamble
%         % The bits are mapped as 1-2*bits since bit 0 is modulated to 1 and bit
%         % 1 is modulated to -1
%         symbol_freq(subcarrier_config~=0) = 1-2*fftshift(bits_preamble);
%         symbol_time = sqrt(num_bins)*ifft(symbol_freq);
%         signal(1+(m-1)*num_bins:m*num_bins) = symbol_time;
%     end
%     signal(1+num_syms_preamble*num_bins:num_syms_preamble*num_bins+cp) = symbol_time(1:cp);
%     % only do one symbol!
% %     symbol_freq(subcarrier_config~=0) = 1-2*fftshift(bits_preamble);
% %     symbol_time = sqrt(num_bins)*ifft(symbol_freq);
% %     signal(1:num_bins) = symbol_time;
% %     signal(1+num_bins:num_bins+cp) = symbol_time(1:cp);
%     % normalize
%     signal = signal./max(abs(signal));
%     rx_signal = rx_signal./max(abs(rx_signal));
%     [c,lags] = xcorr(rx_signal, signal);
%     [t,i] = max(abs(c));
%     disp(lags(i))
%     if t > threshold
%         packet_start_index = lags(i)+1;
%         packet_start_index = packet_start_index + floor((3/4)*cp);
%     else
%         packet_start_index = -1;
%     end
% end
    
%     % adapted from MATLAB OFDMReceiver
%     
%     % generate preamble
%     signal = zeros(1,num_syms_preamble*num_bins);
%     symbol_freq = zeros(1,num_bins);
%     subcarrier_config = ones(1,num_bins);
%     subcarrier_config(convert_bin_index_normal_to_fft(guard_bins,num_bins)) = 0;
%     subcarrier_config(convert_bin_index_normal_to_fft(pilots,num_bins)) = 3; 
%     for m = 1:num_syms_preamble
%         symbol_freq(subcarrier_config~=0) = 1-2*fftshift(bits_preamble);
%         symbol_time = sqrt(num_bins)*ifft(symbol_freq);
%         signal(1+(m-1)*num_bins:m*num_bins) = symbol_time;
%     end
%     %signal(1+num_syms_preamble*num_bins:num_syms_preamble*num_bins+cp) = symbol_time(1:cp);
%     signal = signal./max(abs(signal));
%     rx_signal = rx_signal./max(abs(rx_signal));
%     
%     threshold = 0.6;
%     %L = num_bins;
%     %K = num_bins/4;
%     K = length(signal);
%     %windowLength = ceil(0.5*num_samples + length(signal));
% 
%     % Cross correlate
%     %rWin = rx_signal(1: windowLength-L+K-1);
%     rWin = rx_signal;
%     [Phat,pl] = xcorr(rWin, signal);
%     [Rhat,rl] = xcorr(abs(rWin).^2, ones(K,1)); % Moving average
% 
%     % Remove leading and tail zeros overlaps
%     PhatShort = Phat(ceil(length(Phat)/2):end-K/2+1);
%     RhatShort = Rhat(ceil(length(Rhat)/2):end-K/2+1);
% %     PhatShort = Phat(pl > 0);
% %     RhatShort = Rhat(rl > 0);
% %     PhatShort = PhatShort(1:end-K/2+1);
% %     plot(real(PhatShort))
% %     RhatShort = RhatShort(1:end-K/2+1);
% 
%     % Calculate timing metric
%     M = abs(PhatShort).^2 ./ RhatShort.^2;
%     %plot(M);
% 
%     % Determine start of short preamble. First find peak locations
%     MLocations = find(M > (max(M)*threshold));
% 
%     % Correct estimate to the start of preamble, not its center
%     MLocations = MLocations + K/2;
% 
%     % Pick earliest peak in time
%     packet_start_index = min(MLocations);
%     if isempty(MLocations) || packet_start_index <= 0
%         % No desirable location found
%         packet_start_index = -1; 
%     end


function packet_start_index = packet_detection(rx_signal,P)
    num_bins = P.num_bins;
    num_syms_preamble = P.num_syms_preamble;
    cp = P.cp;
    % DC removal filter
    alpha = 0.975;
    rx_signal = filter([1 -1], [1 -alpha], rx_signal);
    threshold = 6;
    window_size = 4*num_bins;
%     max_ratio = -Inf;
    ratios = zeros(2,numel(rx_signal)-2*window_size+1);
    for j = 1:numel(rx_signal)-2*window_size+1
        ratio = 1/((sum(abs(rx_signal(j:j+window_size-1)).^2))/sum(abs(rx_signal(j+window_size:j+2*window_size-1)).^2)); 
        ratios(1,j)=ratio;
        ratios(2,j)=j+window_size;
%         if ratio > max_ratio
%             packet_start_index = j+window_size;
%             max_ratio = ratio;
%         end
%         ratio_all(j+window_size) = ratio;
    end
    
%     figure
%     hold on
%     plot(abs(rx_signal))
%     plot(0.5*ratio_all/max(ratio_all))
%     plot(packet_start_index,0.4,'*')
%     title('dank')
%     disp(max_ratio)
%     pause(3)
%     if(max_ratio < threshold)
%         packet_start_index = -1;
%     else
%         packet_start_index = packet_start_index + floor((3/4)*cp);
%     end
    [peaks,locs] = findpeaks(ratios(1,:),'MinPeakDistance',num_bins*num_syms_preamble);
    i = find(peaks > threshold, 1);
    if i > 0
        packet_start_index = ratios(2,locs(i))+floor((3/4)*cp);
    else
        packet_start_index = -1;
    end
end
