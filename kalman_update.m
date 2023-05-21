function [X_hat, P_hat] = kalman_update(X_bar, C, R, Y, P_bar)
% Kalman Filter Update Stage
% 
% ==========================Signals=================================
% X_hat, P_hat -> Estimative based on current measurement
% X_bar, P_bar -> Estimative based on previous measurements
% C, R -> State-space model (w/ noise) matrices
% Y -> Measurement (sensor) signal corrupted

% Compute Kalman Gain (teoretical result):
    K = P_bar * (C')*inv(C*P_bar*(C') + R);
% -----------------Matlab warning!!!---------------------------------
% Never use the inverse of a matrix to solve a linear system Ax=b with 
% x=inv(A)*b, because it is slow and inaccurate. Suggested Action:
% Instead of multiplying by the inverse, use matrix right division (/) or 
% matrix left division (\). That is:
% -> Replace inv(A)*b with A\b
% -> Replace b*inv(A) with b/A
%     K = P_bar * (C')/(C*P_bar*(C') + R);
    
% Update:
    X_hat = X_bar + K*(Y - C*X_bar);
    P_hat = P_bar - (K*C*P_bar);
end

