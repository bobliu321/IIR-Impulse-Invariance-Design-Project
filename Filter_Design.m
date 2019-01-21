% Part 5 - Butterworth filter design
% Keith Leung - 301221899
% Bob Liu - 301236133
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delta_l design parameter matrix
delta_l = [0.2929, 0.1, 0.01]
% delta_h design parameter matrix
delta_h = [sqrt(0.1), sqrt(0.1), sqrt(0.01)]


for a = 1:3

% Derivation of the filter order (N) and cutoff frequency (w)
butter_pass = (1 - delta_l(a))^2
butter_stop = (delta_h(a))^2

Inv_butter_pass = (1/butter_pass) - 1
Inv_butter_stop = (1/butter_stop) - 1

butter_filter = Inv_butter_stop / Inv_butter_pass

% The order of filter must be rounded up as an integer
N(a) = ceil((log10(butter_filter) / log10(0.4/0.25)) / 2)

% Frequency in continuous time domain
Omega_c(a) =  (0.25*pi)/(Inv_butter_pass)^(1/(2*N(a)))

% Assuming T = 1, frequency in discrete time domain = frequency in
% continuous time domain
omega_c(a) = Omega_c(a)

% Zeroes, poles and gain of the butterworth transfer function H_c(s) into 
% r, p and k respectively
[r,p,k] = butter(N(a), omega_c(a), 's')

% Forming H_c(s)
s = tf('s')
H_s = 1
for b = 1:numel(p)
    H_s = H_s*(Omega_c(a)/(s - (p(b))))
end

% Coefficients of numerator and denominator polynomials of the butterworth 
% transfer function into x and y
[x,y] = butter(N(a), omega_c(a), 's')

% Partial fraction on x and y matricies and store zeroes, poles and gain 
% into r, p and k respectively
[r,p,k] = residue(x,y)

% Inverse laplace transform of the butterworth transfer function after 
% partial fraction
syms x
syms z
h_t = 0
for d=1:numel(r)
    h_t = h_t + ilaplace((r(d)) / (x - (p(d))))
end

% Z-transform h[n] from inverse laplace transform
H_z = ztrans(h_t)


% Plotting magnitude of DT low pass filter 
%simplify the z-transform
[Numer, Denom] = numden(collect(simplifyFraction(vpa(H_z, 10))))

%obtain numerator and denominator
coeff_Numer = double(coeffs(vpa(Numer)))
coeff_Denom = double(coeffs(vpa(Denom)))
 
figure(a)
freqz(coeff_Numer, coeff_Denom)
title(['Magnitude of Order ', num2str(N(a)), ' Discrete Low Pass Butterworth Filter'])
xlim([0, 1])
ylim([-40, 20])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 6 - Highpass filter design

alpha = -(cos((omega_c(a) + 0.7*pi)/2))/(cos((omega_c(a) - 0.7*pi)/2))
H_zp = subs(H_z, z, -((1 + alpha*(1/z))/(alpha + (1/z))))

end


%Omega vs Theta relation
syms w
figure(4)
f = -atan2((-3/4*sin(w)),(-1-(5/4*cos(w))))
f1 = -atan2(-sin(w),-cos(w))
f2 = -atan2((-3*sin(w)/4),(1-(5*cos(w)/4)))


fplot(f), hold on
fplot(f1), hold on
fplot(f2)
title('Theta vs. Omega')
xlabel('Omega') % x-axis label
ylabel('Theta') % y-axis label
xlim([-pi, pi])
ylim([-pi, pi])