function [r_cfo,r_sfo] = estimate_residual_cfo_sfo(rx_sym,H,P)
    cp = P.cp;
    bits_pilots = P.bits_pilots;
    pilots = P.pilots;
    num_bins = P.num_bins;
    fs = P.fs;
    %
    tx_pilots = 1 - 2*(bits_pilots);
    rx_pilots = rx_sym(convert_bin_index_normal_to_fft(pilots,num_bins))./H(convert_bin_index_normal_to_fft(pilots,num_bins));
    phase_accumulated = angle(rx_pilots./tx_pilots);
    esti_params = regress(phase_accumulated', [ones(1,length(pilots)); pilots]');
    r_cfo = esti_params(1)*fs /(2*pi*(num_bins+cp));
    r_sfo = esti_params(2)*fs /(2*pi*(num_bins+cp)/num_bins);
end
