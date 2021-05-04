function rx_sym_no_channel = correct_channel(rx_sym,H,P)
    guard_bins = P.guard_bins;
    num_bins = P.num_bins;
    %
    rx_sym_no_channel = rx_sym./H;
    rx_sym_no_channel(convert_bin_index_normal_to_fft(guard_bins,num_bins)) = 0;
end
