%% Radar Target Generation and Detection using CFAR

%%

close all

clear all

clc


%% FMCW Radar Specifications 

% Operating frequency of FMCW Radar 
fc = 77e9;

% Speed of light ~ 3e8 m/s
c = physconst('LightSpeed');

% Maximum range = 200 m
Rmax = 200;  

% Range resolution = 1 m
res = 1; 

% Maximum target velocity = 70 m/s = 252 km/h
vmax = 70.0; 

% Velocity resolution = 3 m/s
vres = 3.0;  

%% Target Specifications

% Target's initial position 
rstart = 160.0; 

% Target's velocity (remains constant)
v_t = -20.0; 

%% FMCW Waveform Generation

% Wavelength of Radar carrier signal
lambda = c / fc; 

% Maximum signal round trip time 
rtt = 2 * Rmax / c; 

% Sweep time based on the maximum signal round trip time
Tchirp = 5.5 * range2time(Rmax, c); 

% Bandwidth of the FMCW chirp
B = rangeres2bw(res,c);  

% Slope of the FMCW chirp
slope = B / Tchirp;

% Maximum frequency that the radar has to detect is the sum of the beat frequency corresponding to maximum range and maximum doppler frequency.

% Frequency corresponding to maximum range
f_Rmax = range2beat(Rmax, slope,c); 

% Maximum doppler shift
f_doppler = speed2dop(2*vmax,lambda);

% Maximum beat frequency 
f_max = f_Rmax + f_doppler; 

% Number of chirps in one sequence according to the desired velocity resolution. The whole range from -vmax to +vmax is covered by the doppler cells.
Nd = 2^(nextpow2(2*vmax/vres));

% Number of samples on each chirp
Nr = 2^(nextpow2(Tchirp * f_max) + 1);

% The nyquist criterion needs to be satisfied, so the sampling rate must be more than twice as high as the maximum beat frequency to be resolved.

% Sampling time 
ts = Tchirp / Nr; 

% Sampling Frequency
fs = Nr / Tchirp; 

% Timestamps for running the simulation for every sample on chirp signal sequence
t = linspace(0, Nd * Tchirp, Nr * Nd); 

% Range of the target is updated for each time stamp assuming a constant target velocity
r_t = rstart + v_t * t; 

% Final target range
r_t_final = r_t(length(r_t));

% Time delay between transmitted and received signal for each target range position
td = 2 * r_t / c;  %

% Transmitted signal for each timestamp
Tx = cos(2 * pi * (fc .* t + slope * t.^2 / 2));

% Received signal for each timestamp
Rx = cos(2 * pi * (fc .* (t - td) + slope * (t - td).^2 / 2));

% Beat/mixed signal (found by elementwise matrix multiplication of Tx and Rx signals
Mix = Tx .* Rx;

figure('Name','Transmitted and Received Signals'); 
plot(t, slope.*t, 'm.', t, slope.*(t-td), 'b.');
xlim([0 1e-5]);
legend('Transmitted signal','Received signal');

%% Range Measurement

% Mixed signal vector reshaped into Nr*Nd array
X = reshape(Mix, [Nr, Nd]);

% FFT run on the beat signal and the absolute value of the fft output taken and normalized over the number of data points
X_fft = abs(fft(X, Nr)) / Nr;

% FFT gives a spectrum from -fs/2 to +fs/2 (where fs = sampling frequency) with a symmetry axis going through f = 0. We only use the upper half from f = 0 to +fs/2 and throw away the rest
X_fft = X_fft(1 : Nr / 2 + 1, :);

% Frequency vector from f = 0 to +fs/2 along the chirp or range axis
f_fft = linspace(0, fs / 2, Nr / 2 + 1);

% Range vector
range_fft = c * f_fft * Tchirp / (2 * B);

% Estimated target position (range). There is only one peak.
r_est = range_fft(find(X_fft == max(X_fft), 1, 'first'));

% Plotting the range FFT output
figure('Name','Range measurement by FFT')
title('Fourier Transform of beat signal chirp')
plot(range_fft, X_fft(:, 1))
grid on
set(gca, 'XLim', [0, 200])
xlabel('Range (m)'), ylabel('Signal magnitude')

%% Range Doppler Response

% Range Doppler Map Generation - 2D FFT is applied on the mixed or beat signal to obtain the range doppler response. The first axis corresponds to the length of the chirp signal
% sequence (range axis) and the second axis corrensponds to the number of
% chirp sequences transmitted and received one after another (doppler axis).
% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins.

% Approximate Doppler frequency shift for the target velocity
fdoppler = 2 * v_t / lambda;

% Mix signal vector reshaped into an Nr*Nd arra
Mix = reshape(Mix, [Nr, Nd]);

% 2D FFT run using the FFT size for both dimensions
X_fft2 = fft2(Mix, Nr, Nd);

% Zero doppler frequency component shifted to the center of the spectrum
X_fft2 = fftshift(X_fft2, 2);  

% One side of signal from range dimension taken
X_fft2 = X_fft2(1 : Nr / 2 + 1, 1 : Nd);
RDM = abs(X_fft2);
RDM = 10 * log10(RDM);

% Dopper frequency shift vector along the Doppler axis starting from min(fdoppler) to max(fdoppler)
f_doppler_fft = linspace(-1 / Tchirp / 2, 1 / Tchirp / 2, Nd);

doppler_axis = f_doppler_fft / 2 * lambda;
range_axis = c * f_fft * Tchirp / (2 * B);

% Surface plot used to plot the Range Doppler Response Map
figure('Name', 'Range Doppler Map')
surf(doppler_axis, range_axis, RDM, 'LineStyle', ':', 'LineWidth', 0.5);
xlabel('Doppler or axis (m/s)')
ylabel('Range axis (m)')
zlabel('Signal level (dB)')
title('Range Doppler Map')
colormap("turbo")
set(gcf, 'Color', 'w', 'Position', [676, 549, 720, 480])


%% 2D CFAR implementation

% 2D CFAR implemented by sliding a window of training and guard cells around
% a cell under test (CUT) on the Range Doppler Map RDM[x,y] obtained from
% the 2D FFT.

% Number of Training Cells selected in both the dimensions
Tr = 16;  % range dimension
Td = 8;  % doppler dimension

% Number of Guard Cells selected in both dimensions
Gr = 8;  % range dimension
Gd = 4;  % doppler dimension

% Threshold offset by SNR value in dB
SNR_offset_dB = 8;

% Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);

% Design a loop such that slides the cell under test (CUT) across the Range
% Doppler Map by giving margins at the edges for Training and Guard Cells.
% For every iteration sum the signal level within all the training cells.
% To sum up the signal values convert the value from logarithmic to linear
% using db2pow function. Average the summed values over all the training
% cells used. After averaging convert it back to logarithimic using pow2db.
% Further add the offset to it to determine the threshold. Next, compare 
% the signal level in CUT with this threshold. If the CUT signal level >
% threshold assign the CUT a value of 1 for a detection, else set it to 0
% for no detection.

% The process above will generate a thresholded block, which is smaller 
% than the Range Doppler Map as the CUT cannot be located at the edges of
% matrix. Hence, few cells will not be thresholded. To keep the map size 
% same set those values to 0.
cfar_2d = zeros(size(RDM));

% Width of the 2D CFAR sliding window in Range and Doppler dimension
wr = 2 * (Gr + Tr) + 1;
wd = 2 * (Gd + Td) + 1;

% 2D array to hold the threshold values
threshold_cfar = zeros(Nr / 2 + 1 - 2 * (Tr + Gr), Nd - 2 * (Td + Gd));

% 2D array to hold the final signal after thresholding
sig_cfar2D = zeros(Nr / 2 + 1 - 2 * (Tr + Gr), Nd - 2 * (Td + Gd));

% Generate 2D mesh grid the cfar threshold and filtered signal
[X_cfar,Y_cfar] = meshgrid((Td + Gd) : 1 : (Nd - (Td + Gd) - 1), ...
    (Tr + Gr) : 1 : (Nr / 2 + 1 - (Tr + Gr) - 1));

% Window slid across the rows of the 2D FFT RDM array where (i, j) 
% is the lower left starting point of the 2D sliding window
for i = 1 : (Nr/2+1 - wr + 1)
    
    % Window slid across the columns of the 2D FFT RDM array where (i, j) 
    % is the lower left starting point of the 2D sliding window.
    for j = 1 : (Nd - wd + 1)
        
        % Noise threshold determined by measuring the noise level in the
        % training cells of the sliding window within the 2D FFT
        % converted from logarithmic to linear signal power.
        noise_level = ...
            sum(sum(db2pow(RDM(i : i + wr - 1, j : j + wd - 1)))) - ...
            sum(sum(db2pow(RDM(i + Tr : i + Tr + 2 * Gr + 1, ...
            j + Td : j + Td + 2 * Gd + 1))));
        
        % Number of training cells
        NT = wr * wd - (2 * Gr + 1) * (2 * Gd + 1);
        
        % Noise threshold determined by taking the average of summed noise
        % over all training cells, and converting it back to logarithmic signal
        % values and adding the logarithmic SNR offset.
        threshold = pow2db(noise_level / NT) + SNR_offset_dB;
        threshold_cfar(i, j) = threshold;
          
        % Signal within the CUT measured
        signal = RDM(Tr + Gr + i, Td + Gd + j);
        
        if (signal < threshold)
            sig_cfar2D(i, j) = 0;
            cfar_2d(Tr + Gr + i, Td + Gd + j) = 0;
        else
            sig_cfar2D(i, j) = signal;
            cfar_2d(Tr + Gr + i, Td + Gd + j) = 1;
        end        
        
    end
    
end

% Display the CFAR2D output using the surf() function like we did for Range
% Doppler Response output.
figure('Name', '2D CFAR on Range Doppler Map')
surf(doppler_axis, range_axis, cfar_2d, 'LineStyle', ':', 'LineWidth', 0.5)
colorbar
xlabel('Doppler or velocity axis [m/s]')
ylabel('Range axis [m]')
zlabel('Signal level [dB]')
title('2D CFAR Output on Range Doppler Map')
set(gcf, 'Color', 'w', 'Position', [676, 549, 720, 480])
 
% % Find the estimated target range and velocity in 2D CFAR map (assumption:
% % only one target peak exists in this simulation).
% [rows, cols, vals] = find(cfar2D == 1);
% r_est = range_axis(round((min(rows) + max(rows)) / 2));
% v_t_est = doppler_axis(round((min(cols) + max(cols)) / 2));
% 
% % Print true and estimated target range and velocity
% fprintf('\nEstimated target range and velocity from 2D CFAR:\n');
% fprintf('-------------------------------------------------------------\n');
% fprintf('True average target range   r_t_mean = %6.4e m\n', ...
%     (rstart + r_t_final) / 2);
% fprintf('Estimated target range       r_t_est = %6.4e m\n', r_est);
% fprintf('\n');
% fprintf('True target velocity             v_t = %6.4e m/s\n', v_t);
% fprintf('Estimated target velocity    v_t_est = %6.4e m/s\n', v_t_est);
% fprintf('-------------------------------------------------------------\n');
% fprintf('\n');

