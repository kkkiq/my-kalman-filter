function [X_bar, P_bar] = kalman_propagation(X_hat, A, B, U, G, Q, Ts, P_hat)
% Kalman Filter Propagation Stage using Euler1 numerical integration
% 
% ==========================Signals=================================
% X_bar, P_bar -> Estimative based on previous measurements
% X_hat, P_hat -> Estimative based on current measurement (delayed by Z⁻¹)
% A, B, G, Q -> State-space model (w/ noise) matrices
% U -> Input signal (ideal)
% Ts -> Sampling period

    X_bar = (A*X_hat + B*U)*Ts + X_hat;
    P_bar = (A*P_hat + P_hat*(A') + G*Q*(G'))*Ts + P_hat;
end

