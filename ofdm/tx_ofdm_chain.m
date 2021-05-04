function [signal] = tx_ofdm_chain(bits_data,P)
    % parameter unpack
    num_samples = P.num_samples;
    num_bins = P.num_bins;
    guard_bins = P.guard_bins;
    pilots = P.pilots;
    num_syms_preamble = P.num_syms_preamble;
    bits_preamble = P.bits_preamble;
    cp = P.cp;
    num_syms_data = P.num_syms_data;
    num_bins_data = P.num_bins_data;
    bits_pilots = P.bits_pilots;
    % tx chain
    signal = zeros(1,num_samples);
    symbol_freq = zeros(1,num_bins);
    subcarrier_config = ones(1,num_bins);
    subcarrier_config(convert_bin_index_normal_to_fft(guard_bins,num_bins)) = 0;
    subcarrier_config(convert_bin_index_normal_to_fft(pilots,num_bins)) = 3; 
    % preamble symbols
    for m = 1:1:num_syms_preamble
        symbol_freq(subcarrier_config~=0) = 1-2*fftshift(bits_preamble);
        symbol_time = sqrt(num_bins)*ifft(symbol_freq);
        signal(1+(m-1)*num_bins:m*num_bins) = symbol_time;
    end
    signal(1+num_syms_preamble*num_bins:num_syms_preamble*num_bins+cp) = symbol_time(1:cp);
    % data symbols
    for m = 1:1:num_syms_data
        symbol_freq(subcarrier_config==1) = 1-2*fftshift(bits_data(1+(m-1)*num_bins_data:m*num_bins_data));
        symbol_freq(subcarrier_config==3) = 1-2*fftshift(bits_pilots);
        symbol_time = sqrt(num_bins)*ifft(symbol_freq);
        signal(1+num_syms_preamble*num_bins+cp+(m-1)*(num_bins+cp):m*(num_bins+cp)+num_syms_preamble*num_bins+cp) = [symbol_time, symbol_time(1:cp)];
    end
end
