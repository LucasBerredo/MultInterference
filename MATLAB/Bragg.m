%% Script principal
lam = 1;          % Central wavelength
n_High = 4;       % High index
n_Low = 1.5;      % Low index
n_inc = 1;        % Incident medium
n_sub = 1.5;      % Substrate
m = 20;           % Number of pairs
theta_i0 = 0;     % Angle of incidence (rad)
mesh_len = 1000;  % Mesh length

bragg(lam,n_High,n_Low,n_inc,n_sub,m,theta_i0,mesh_len)
%%
function bragg(lam, n_High, n_Low, n_inc, n_sub, m, theta_i0s, mesh_len)

%% Setup
n = [n_inc, repmat([n_High n_Low],1,m), n_sub];     % Indices
d = repmat([lam/(4*n_High) lam/(4*n_Low)],1,m);     % lambda/4 thickness
mesh = linspace(0.5*lam,1.5*lam,mesh_len);         % Mesh

angles = theta_i0s(:);   % Ensure column vector

figure('Position',[100 100 800 500])
hold on

for k = 1:length(angles)

    theta = angles(k);

    [~,~,~,~,~,r,~] = MultInterference(1, theta, n, d, mesh, false);
    R = abs(r).^2;

    %% Theoretical values
    BW = 4*asin((n_High-n_Low)/(n_High+n_Low))/pi;
    lam_BW_min = lam/(1+BW/2);
    lam_BW_max = lam/(1-BW/2);

    ratio_idx = n_sub*((n_High/n_Low)^(2*m));
    R_max = ((n_inc-ratio_idx)/(n_inc+ratio_idx))^2;

    %% Plotting
    figure('Position',[100 100 800 500])

    plot(mesh,R,'b','LineWidth',2)
    hold on

    yline(R_max,'r--','Theoretical R_{max}')
    xline(1.0,'g:','LineWidth',3,'Central Wavelength (\lambda_0)')

    ylim([0 1.1])

    xlabel('Normalized Wavelength (\lambda / \lambda_0)')
    ylabel('Reflectance')
    title('Bragg Mirror Reflectance Spectrum')

    legend
    grid on

end

end




%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
