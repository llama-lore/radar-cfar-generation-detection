# Generation and CFAR Detection of Target Using MATLAB


A Radar Engineering Project

## Project Summary

Constant False Alarm Rate (CFAR) is an adaptive threshold whose level is determined by the clutter and / or noise in the vicinity of the radar echo. In radar and other sensing applications, the ability to detect weak targets amidst strong interference is paramount. CFAR achieves this by adaptively adjusting its detection threshold based on the local statistical properties of the background clutter, ensuring a constant probability of false alarm across different clutter environments. The objective of this project is to explore how CFAR can be effectively implemented using MATLAB for both the generation of target signals and the subsequent detection of these targets in noisy and cluttered environments.


<p align="center"><img width = "700" src="https://github.com/user-attachments/assets/943a039f-213e-470f-91ed-bc05a4123994"</p>

The simulation involves the implementation of the following steps on MATLAB-
1.	First, the specifications of the FMCW radar are defined. This includes the operating frequency (fc), the speed of light (c), maximum range (Rmax), range resolution (res), assumed maximum target velocity (vmax) and the velocity resolution (vres).
2.	The target’s specifications i.e. the target’s velocity (v_t) and initial position (rstart) are then defined. They are verified through different methods throughout the simulation. 
3.	The chirp time (Tchirp), the slope of the chirp (Slope), the bandwidth (B), the beat frequency corresponding to maximum range (beat_Rmax) and the doppler shift (beat_doppler) are calculated. This helps calculate the maximum beat frequency (beat_max) [beat_max = beat_Rmax + beat_doppler].
4.	 The number of chirps in one sequence and the number of samples on each chirp are defined to find the sampling frequency. The sampling rate must be at least twice as much as the maximum desired beat frequency.
5.	Timestamps are defined for every sample on the chirp sequence and the range of the target is updated for each time stamp assuming constant target velocity.
6.	Transmitted and Received signals are defined for each timestamp.

<p align="center"><img width = "500" src="https://github.com/user-attachments/assets/61c0f108-23b2-40af-95ea-90b8f6de8e5e"</p>

7.	With the waveform parameters, the target movement is simulated by calculating the beat signal (mixed signal) for every timestamp. The beat signal gives the correct range (approximate) by using FFT.

<p align="center"><img width = "500" src="https://github.com/user-attachments/assets/14695163-e433-4b46-aec3-16cf3a752c86"</p>


8.	The FFT is implemented on the mixed signal, and the result is plotted. Correct implementation generates a single peak at the correct target position (range). 
9.	The range doppler response is then found by using 2D FFT on the chirp signal. The first axis corresponds to the length of the chirp signal (range axis). The second axis corresponds to the number of chirp sequences that are transmitted and received (doppler axis). The output is an image that has a response in the range and doppler FFT bins.
10.	The output of the 2D FFT is observed using a surface plot. It is verified that the peak corresponds to the range originally considered.

<p align="center"><img width = "500" src="https://github.com/user-attachments/assets/d0e25831-fb06-46c7-b5b2-2b806a7f1bf1"</p>

11.	The last step is the implementation of the 2D CFAR process on the range doppler map. This can suppress the noise and separate the target signals. This is done by defining training cells and guard cells. The output of the CFAR process indicates a range and velocity estimate of the target.


