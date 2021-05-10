function [mag, fftCoeff, freqList] = plot_fft(data, sampleRate)
    data = data(:);
    fftCoeff = fftshift(fft(data), 1);
    mag = abs(fftCoeff);
    N = length(data);
    % calculate frequency spacing
    df = sampleRate / N;
    % calculate unshifted frequency vector
    f = (0:(N-1))*df;
    % move all frequencies that are greater than fs/2 to the negative side of the axis
    f(f >= sampleRate/2) = f(f >= sampleRate/2) - sampleRate;
    % freq are aligned with one another; if you want frequencies in strictly
    % increasing order, fftshift() them
    freqList = fftshift(f);
    if(nargout == 0)
        figure;
        plot(freqList/1e6, 20*log(mag),'k');
        xlabel('Frequency (MHz)');
        ylabel('Magnitude (dB)');
        %set(gca, 'FontSize', 15);
    end
end
