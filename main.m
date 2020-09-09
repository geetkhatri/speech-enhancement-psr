% Description: 
% Speech enhancement using Wiener filtering and phase reconstruction of 
% the pitch-synchronous representation (PSR) of the speech signal
%
% Usage:
% Run the script. Provide the path to the audio file when prompted.
% The samples of the enhanced speech signal will be stored as y.
%
% Author:
% Geet Khatri
% https://github.com/geetkhatri/

% prompt the user to enter filename
filename = input('\n\nFilename: ', 's');
while exist(filename, 'file') ~= 2
    filename = input('File does not exist. Try again: ', 's');
end

info = audioinfo(filename);
if info.NumChannels ~= 1
    error('Audio must be mono.')
end

tic         % for timing

% Input signal
[x, fs] = audioread(filename);

fsr = 16000;                        % frequency of resampled signal
M_default = 0.032 * fsr;            % min. window length
overlap = 0.875;                    % overlap between adjacent frames
H = (1 - overlap) * M_default;      % hop size
K = 6;                              % ratio M/(fs*T0)

% Resampling
if fs ~= fsr
    [D, nframes] = rat(fs / fsr);   % decimation and interpolation factors
    x = resample(x, nframes, D);    % resampled signal
    fs = fsr;
end

% Add white Gaussian noise to signal
SNR = 30;                   % signal-to-noise ratio
x = awgn(x, SNR);           % noisy signal

% Fundamental frequency detection
p.fres = 1.81 * fs / M_default;       % frequency resolution (for fxpefac)
tinc = H/fs;                   % time increment between consecutive frames
cd voicebox
[f0, tv, pv, ~] = fxpefac(x, fs, tinc, 'G', p);  % fundamental frequencies
nframes = length(tv);                % number of frames
voiced = zeros(nframes, 1);          % tells whether a frame is voiced
for i = 1 : nframes
    if pv(i) >= 0.9
        voiced(i) = 1;
    end
end
cd ../

% Wiener filtering
IS = M_default;     % initial silence
i = 2;
while voiced(i) ~= 1
    IS = IS + H;
    i = i + 1;
end
x_noisy = x;        % save noisy signal for comparison
[esTSNR, x] = WienerNoiseReduction(x, fs, IS);  % filtered signal

% Computation of window sizes
min_overlap = 0.73;              % minimum overlap for Hamming window
M_min = round(H / (1 - min_overlap));           % minimum window size
M_max = 1024;                                   % maximum window size
M = zeros(1, nframes);           % array of window sizes
for l = 1 : nframes
    if voiced(l) == 1
        M(l) = round(K * fs / f0(l));
        if M(l) < M_min
            M(l) = M_min;
        elseif M(l) > M_max
            M(l) = M_max;
        end
    else
        M(l) = M_default;
    end
end

% Computation of STFT
X = cell(1, nframes);         % STFT cell array
Mh = zeros(1, nframes);       % half window size
w = cell(1, nframes);         % window cell array
offx = M_default / 2;         % offset in x
for i = 1 : nframes
    % extract frame of input data
    if mod(M(i), 2) == 0    % even M
        Mh(i) = M(i) / 2;
        xf = x(offx - Mh(i) + 1 : offx + Mh(i));
    else
        Mh(i) = (M(i) - 1) / 2;     % odd M
        xf = x(offx - Mh(i) : offx + Mh(i));
    end
    w{i} = hamming(M(i), 'periodic');   % window function
    xfw = xf .* w{i};             % apply long window to current frame
    X{i} = fft(xfw);         % STFT for frame l
    offx = offx + H;            % advance pointer by hop-size H
end

% Phase reconstruction
Y = X;                    % initialise spectrum of enhanced signal
for l = 2 : nframes
    
    % voiced sound
    if voiced(l) == 1
       
       magY = abs(X{l});
       phaY = princ(angle(X{l}));
       if mod(M(l), 2) == 0
           nbands = 1 + M(l) / 2;
       else
           nbands = (1 + M(l)) / 2;
       end
       phaY = phaY(1 : nbands);
        
       % number of harmonics
       nharm = floor(fs ./ (2 * f0(l)));
       
       % harmonic indices
       r = (1 : nharm).';
       
       % harmonic frequencies
       fh = r * f0(l);
       
       % bins that contain harmonics
       bh = r * K;

       % not onset of voiced sound
       if voiced(l - 1) ~= 0
           % normalised fk in radians
           wh = 2 * pi * fh / fs;
       
           % phase reconstruction along time
           if nharm_prev >= nharm
               phaY(bh + 1) = princ(phaY_prev(bh + 1) + wh * H);
           else
               r_prev = (1 : nharm_prev).';
               bh_prev = r_prev * K;
               wh_prev = 2 * pi * r_prev * f0(l - 1) / fs;
               phaY(bh_prev + 1) = princ(phaY_prev(bh_prev + 1) ...
                                                    + wh_prev * H);
           end
       end
       
       % extent of leakage of harmonic components into neighboring bands
       delk = ceil(K / 2);
       
       % phase reconstruction along frequency
       dtft_phase = angle(fft(w{l}));
       for i = [-delk : -1, 1 : delk]
           diff = sign(i) * dtft_phase(abs(i) + 1);
           for j = 1 : nharm
               if 1 + bh(j) + i <= nbands
                   phaY(1 + bh(j) + i) = princ(phaY(1 + bh(j)) + diff);
               end
           end
       end
       
       % save for next iteration
       nharm_prev = nharm;
       phaY_prev = phaY;
       
       % reconstruct phase spectrum for remaining frequencies
       phaY2 = zeros(M(l), 1);
       phaY2(1 : nbands) = phaY;
       temp = -flip(phaY);    
       if mod(M(l), 2) == 0
           temp = temp(2 : end - 1);
       else
           temp = temp(1 : end - 1);
       end
       phaY2(nbands + 1 : end) = temp;
       Y{l} = magY .* exp(1j * phaY2);
       
    end
    
end

% Reconstruction of signal using overlap-add method
leny = max(length(x), M_default / 2 + (nframes - 1) * H + Mh(end));
                                          % length of reconstructed signal
y = zeros(leny, 1);             % reconstructed signal
NF = zeros(leny, 1);            % normalization function
offy = M_default / 2;           % offset in y
for i = 1 : nframes
    if mod(M(i), 2) == 0
        samples = offy - Mh(i) + 1 : offy + Mh(i);
    else
        samples = offy - Mh(i) : offy + Mh(i);
    end
    yfw = ifft(Y{i}) .* w{i};
    y(samples) = y(samples) + yfw;    
    NF(samples) = NF(samples) + (w{i} .^ 2);
    offy = offy + H;
end
y = y ./ NF;                % normalise y
y = real(y);                % ignore imaginary parts

toc         % for timing

% N = 1024;
% w1 = hamming(M_default);
% 
% [X1, f1, t1] = spectrogram(x, w1, M_default - H, N, fs);
% figure
% surf(t1, f1 / 1000, princ(angle(X1)), 'EdgeColor', 'None')
% view(2)
% title('Phase of noisy signal')
% xlabel('Time (secs)')
% ylabel('Frequency (kHz)')
% 
% [Y1, f1, t1] = spectrogram(y, w1, M_default - H, N, fs);
% figure
% surf(t1, f1 / 1000, princ(angle(Y1)), 'EdgeColor', 'None')
% view(2)
% title('Phase of reconstructed signal')
% xlabel('Time (secs)')
% ylabel('Frequency (kHz)')