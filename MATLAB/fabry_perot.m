%% Script principal
lam0 = 1.0;
n_High = 3.5;
n_Low = 1.5;
n_cavity = 3.5;
n_inc = 1.0;
n_sub = 1.5;
m = 4;
cav_mult = 2.0;
theta_i0 = 0.0;
mesh_len = 2000;

fabry_perot(lam0,n_High,n_Low,n_cavity,n_inc,n_sub,m,cav_mult,theta_i0,mesh_len)
%%

function fabry_perot(lam0, n_High, n_Low, n_cavity, n_inc, n_sub, m, cav_mult, theta_i0s, mesh_len)

%% Mirror 1
n_mirror = repmat([n_High n_Low],1,m);
d_mirror = repmat([lam0/(4*n_High) lam0/(4*n_Low)],1,m);

%% Mirror 2 (reversed)
n_mirror2 = fliplr(n_mirror);
d_mirror2 = fliplr(d_mirror);

%% Stack
n = [n_inc n_mirror n_cavity n_mirror2 n_sub];
d = [d_mirror cav_mult*lam0/(4*n_cavity) d_mirror2];

mesh = linspace(0.7*lam0,1.3*lam0,mesh_len);

%% Interference calculation
[~,~,Pi,~,Pt,r,~] = MultInterference(1.0, theta_i0s, n, d, mesh, false);

R = abs(r).^2;
T = Pt ./ Pi;

%% Plotting
figure('Position',[100 100 800 500])
hold on

plot(mesh/lam0, T,'b','LineWidth',2)
plot(mesh/lam0, R,'r--','LineWidth',1.5)

xline(1.0,'g:','First Order Resonance (\lambda_1)','LineWidth',1.5)

ylim([-0.05 1.05])
xlim([0.7 1.3])

xlabel('Normalized Wavelength (\lambda / \lambda_0)')
ylabel('Power Fraction')

title('Fabry-Pérot Resonant Cavity Spectrum')

legend('Transmittance (T)','Reflectance (R)','Location','east')
grid on

end




%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
