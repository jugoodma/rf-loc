function [bits_data] = rx_ofdm_chain(rx_signal,P)
    if ~isrow(rx_signal)
        rx_signal = transpose(rx_signal);
    end
    % parameter unpack
    cp = P.cp;
    pilots = P.pilots;
    guard_bins = P.guard_bins;
    num_bins = P.num_bins;
    num_bins_data = P.num_bins_data;
    num_syms_data = P.num_syms_data;
    %
    % bits_data = [];
    bits_data = zeros(1,num_syms_data*num_bins_data);
    subcarrier_config = ones(1,num_bins);
    subcarrier_config(convert_bin_index_normal_to_fft(guard_bins,num_bins)) = 0;
    subcarrier_config(convert_bin_index_normal_to_fft(pilots,num_bins)) = 3; 
    current_index = 1;
    % packet detection
    packet_start = packet_detection(rx_signal(current_index:end),P);
    current_index = packet_start + 4*num_bins;
    % cfo
    cfo = estimate_cfo(...
        rx_signal(current_index:current_index+num_bins-1),...
        rx_signal(current_index+num_bins:current_index+num_bins+num_bins-1),...
        P...
    );
    rx_signal(current_index:end) = correct_cfo(rx_signal(current_index:end),cfo,P);
    current_index = current_index + 2*num_bins;
    % channel
    H = estimate_channel(...
        rx_signal(current_index:current_index+num_bins-1),...
        rx_signal(current_index+num_bins:current_index+2*num_bins-1),...
        P...
    );
    current_index = current_index + 2*num_bins+cp;
    % decoding
    for m = 1:1:num_syms_data
        symt = rx_signal(current_index:current_index+num_bins-1);
        symf = (1/sqrt(num_bins))*fft(symt);
        [r_cfo,r_sfo] = estimate_residual_cfo_sfo(symf,H,P);
        H = correct_residual_cfo_sfo(H,r_cfo,r_sfo,P);
        symf = correct_channel(symf,H,P);
        bits = symf(subcarrier_config==1) < 0;
        % bits_data = [bits_data, fftshift(bits)];
        bits_data(1+(m-1)*num_bins_data:m*num_bins_data) = fftshift(bits);
        current_index = current_index + num_bins + cp;
    end
end
