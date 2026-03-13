function [Hr, Ht, Pi, Pr, Pt, R, T] = MultInterference(Hi, theta_i0, n, d, lam, ind)

%% Check whether to vectorize
lam_arr = lam(:);              % Ensure column vector
num_lam = length(lam_arr);

%% Variable initializations
theta = zeros(size(n));
theta(1) = theta_i0;

%% Finding angles
for i = 2:length(n)
    theta(i) = asin(n(i-1) * sin(theta(i-1)) / n(i));
end

%% Finding kz
kz = 2*pi * (n .* cos(theta));
kz = kz ./ lam_arr;   % implicit expansion (MATLAB >= R2016b)

%% Loop to find M
M = repmat(eye(2),1,1,num_lam);

for i = 2:length(n)

    n_cos_prev = n(i-1) * cos(theta(i-1));
    n_cos_curr = n(i)   * cos(theta(i));

    t = 2 * n_cos_prev / (n_cos_prev + n_cos_curr);
    r = (n_cos_prev - n_cos_curr) / (n_cos_prev + n_cos_curr);

    D = (1/t) * [1 r; r 1];

    for k = 1:num_lam
        M(:,:,k) = M(:,:,k) * D;
    end

    if i < length(n)

        phi = kz(:,i) * d(i-1);

        for k = 1:num_lam
            P = [exp(1i*phi(k)) 0; 0 exp(-1i*phi(k))];
            M(:,:,k) = M(:,:,k) * P;
        end

    end
end

%% Final outputs
T = zeros(num_lam,1);
R = zeros(num_lam,1);

for k = 1:num_lam
    T(k) = 1 / M(1,1,k);
    R(k) = M(2,1,k) * T(k);
end

Hr = R * Hi;
Ht = T * Hi;

Pi_scalar = n(1)*real(cos(theta(1)))*abs(Hi)^2;
Pi = Pi_scalar * ones(num_lam,1);

Pr = n(1)  * real(cos(theta(1)))  * abs(Hr).^2;
Pt = n(end)* real(cos(theta(end))) * abs(Ht).^2;

if ind
    fprintf('Hr = %.4f * exp(i*%.4f)\n', abs(Hr(1)), angle(Hr(1)));
    fprintf('Ht = %.4f * exp(i*%.4f)\n', abs(Ht(1)), angle(Ht(1)));
end

if isscalar(lam)
    Hr = Hr(1);
    Ht = Ht(1);
    Pi = Pi(1);
    Pr = Pr(1);
    Pt = Pt(1);
    T  = T(1);
end

end

%[appendix]{"version":"1.0"}
%---
