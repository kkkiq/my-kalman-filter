# My Kalman Filter
A simple Kalman Filter example made in MATLAB.

Task given at Signals and Systems class (2022.1) at Instituto Federal Fluminense (IFF). It uses a linearized PMSM electric motor model as plant, and its goal is to estimate the current states after the Clark and Park transforms (d-q frame) together with its angular velocity, considering that the model is in presence of noises and disturbances with known parameters (mean and variance).

To run the simulation, download the code and open the folder in MATLAB. Then, run "PMSM_Kalman_Filter.m". Its stages (propagation and update) are coded in separate functions, called by a custom block inside the SIMULINK simulation.

![gaussian_estimation_example](https://github.com/kkkiq/my-kalman-filter/assets/85909385/1825162c-e7b3-4e76-bc7d-41bc1fcac4da)

Estimation provided by the Kalman filter. It's height consists of the probability density function of each sample.
