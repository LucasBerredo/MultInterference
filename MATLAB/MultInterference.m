function [Hr,Ht,Pi,Pr,Pt,R,T] = MultInterference(Hi, theta_i0, n, d, lambda, ind)

% ARGUMENTOS DE ENTRADA
% Hi        <- amplitud compleja del haz incidente
% theta_i0  <- ángulo de la primera incidencia en grados
% n         <- vector de N componentes con los índices de refracción (incidente, capas, sustrato)
% d         <- vector de N-2 componentes con los ESPESORES de las capas intermedias en metros
% lambda    <- longitud de onda del haz incidente en metros
% ind       <- indicador de display booleano
%
% Nota: Este código asume polarización TE (Transverse Electric / s-polarization)

% Convertimos grados a radianes para soportar números complejos en TIR (Total Internal Reflection)
theta_i0_rad = theta_i0 * pi / 180;

% Inicializaciones
theta = zeros(1, length(n));
theta(1) = theta_i0_rad;
kz = zeros(1, length(n));
kz(1) = 2*pi*cos(theta(1))*n(1)/lambda;

% Pre-calculamos todos los ángulos y vectores de onda kz
for j = 1:length(n)-1
    theta(j+1) = asin( n(j)*sin(theta(j)) / n(j+1) );
    kz(j+1) = 2*pi*n(j+1)*cos(theta(j+1)) / lambda;
end

M = eye(2);

% Bucle para calcular la matriz total M (Dynamical * Propagation)
for j = 1:length(n)-1
    
    % Matriz Dinámica (Interfaz j a j+1)
    t = 2*n(j)*cos(theta(j)) / ( n(j)*cos(theta(j)) + n(j+1)*cos(theta(j+1)) );
    r = ( n(j)*cos(theta(j)) - n(j+1)*cos(theta(j+1)) ) / ( n(j)*cos(theta(j)) + n(j+1)*cos(theta(j+1)) );
    
    D = (1/t) * [1, r; 
                 r, 1];
                 
    M = M * D;
    
    % Matriz de Propagación (Solo si no es el sustrato final)
    if j < length(n) - 1
        % d(j) es el espesor de la capa intermedia j
        phi = kz(j+1) * d(j);
        P = [exp(1i*phi), 0; 
             0, exp(-1i*phi)];
             
        M = M * P;
    end
end

% Coeficientes de amplitud de reflexión y trasmisión
T = 1 / M(1,1);
R = M(2,1)*T;

% Amplitudes de campo reflejadas y trasmitidas
Hr = R*Hi;
Ht = T*Hi;

% Potencias
Pi = n(1)   * real(cos(theta(1)))   * abs(Hi)^2;
Pr = n(1)   * real(cos(theta(1)))   * abs(Hr)^2;
Pt = n(end) * real(cos(theta(end))) * abs(Ht)^2;

if ind
    % Display en forma exponencial compleja
    disp(['Hr = ', num2str(abs(Hr)), ' * exp(i*', num2str(angle(Hr)), ')']);
    disp(['Ht = ', num2str(abs(Ht)), ' * exp(i*', num2str(angle(Ht)), ')']);
end

end

%[appendix]{"version":"1.0"}
%---
