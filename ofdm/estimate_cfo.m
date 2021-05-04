function cfo = estimate_cfo(rx_sym1,rx_sym2,P)
    fs = P.fs;
    cfo = -(fs/(2*pi*length(rx_sym1)))*angle(sum(conj(rx_sym2).*rx_sym1));
end
