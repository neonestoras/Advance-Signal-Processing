%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD estimate i.e. Periodogram   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For Section 3
% function named pgm that 
% calculates the periodogram (output)
% of a sequence x = [x[1],x[2],...,x[N]] (input)

function [periodogram, normalised_frequency] = pgm(x)
    N = length(x); %Length of input sequence
    normalised_frequency = 0:1/N:(N-1)/N;
 
    periodogram = zeros(N,1); %Initialise periodogram vector
    fft_inputseq=[]; n=[0:1:N-1]';
    fft_inputseq = sum(exp(-1i*2*pi*n*normalised_frequency)*x',2); %Fast Fourier Transform of input sequence,x
    periodogram =(( abs(fft_inputseq).^2 ).*(1/N))'; %obtaining periodogram, as in given equation 
end