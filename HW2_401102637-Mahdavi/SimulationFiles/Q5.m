%% PID Synthesis through Transfer Function Matching - Verification
clc; clear;
s = tf('s');

% Original controller (b)
n2=5; n1=12; n0=8; d2=2; d1=7;
C_orig = (n2*s^2 + n1*s + n0) / (d2*s^2 + d1*s);

% Compute PID parameters using derived formulas
Kp = n1/d1 - n0*d2/d1^2;
Ki = n0/d1;
Kd = n2/d1 - n1*d2/d1^2 + n0*d2^2/d1^3;
tauD = d2/d1;

fprintf('Kp = %s = %.4f\n', rats(Kp), Kp);
fprintf('Ki = %s = %.4f\n', rats(Ki), Ki);
fprintf('Kd = %s = %.4f\n', rats(Kd), Kd);
fprintf('tauD = %s = %.4f\n', rats(tauD), tauD);

% Build standard PID
C_PID = Kp + Ki/s + Kd*s/(tauD*s+1);

% Compare step responses (just a sanity check)
t = 0:0.01:5;
[y1,~] = step(C_orig, t);
[y2,~] = step(C_PID, t);
err = norm(y1-y2);
fprintf('Max difference between responses: %e\n', err);

% Check equivalence of transfer functions
diff_tf = minreal(C_orig - C_PID);
[~, den] = tfdata(diff_tf);
if den{1} == 0
    disp('The transfer functions are identical.');
else
    disp('Difference exists.');
end