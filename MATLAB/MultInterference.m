function [Hr,Ht,Pi,Pr,Pt] = MultInterference(Hi,theta_i0, n, d, lambda,ind)

%ARGUMENTOS DE ENTRADA
% Hi        <- amplitud compleja del haz incidente
% theta_i0  <- ángulo de la primera incidencia
% n         <- vector de N+1 componentes con los índices de refracción
% d         <- vector de N componentes con las N distancias entre capas
% lambda    <- longitud de onda del haz incidente
%
%ARGUMENTOS DE SALIDA
% Hr        <- haz reflejado
% Ht        <- haz trasmitido
% Pi        <- potencia del haz incidente
% Pr        <- potencia del haz reflejado
% Pt        <- potencia del haz trasmitido
% ind       <- indicador de display booleano: 1 muestra el campo como exponenciales imaginarias; 0 no lo muestra
%
%Inicializaciones
theta = zeros(1,length(n));
theta(1) = theta_i0;
kz = zeros(1,length(n));
kz(1) = 2*pi*cos(theta(1))*n(1)/lambda;
M = eye(2);

%Bucle para calcular la matriz total M
for j=1:length(n)-1
theta(j+1) = asin ( n(j)*sin(theta(j)) / n(j+1) );
kz(j+1) = 2*pi* n(j+1)* cos(theta(j+1)) / lambda;

t = 2*n(j)*cos(theta(j)) / ( n(j)*cos(theta(j)) + n(j+1)*cos(theta(j+1)) );
r = t * ( kz(j)-kz(j+1) ) / ( 2*kz(j) );

M = M * (1/t)*[
    exp(1i*(kz(j+1)-kz(j))*d(j)),        r*exp(-1i*(kz(j+1)-kz(j))*d(j));

    r*exp(1i*(kz(j+1)-kz(j))*d(j)),      exp(-1i*(kz(j+1)-kz(j))*d(j))];
end

%Coeficientes de reflexión y trasmisión
t = 1 / M(1,1);
r = M(2,1)*t;
%Amplitudes de campo reflejadas y trasmitidas
Hr = r*Hi;
Ht = t*Hi;
%Potencias
Pi = n(1)*      cos(theta(1))*    abs(Hi)^2;
Pr = n(1)*      cos(theta(1))*    abs(Hr)^2;
Pt = n(end)*    cos(theta(end))*  abs(Ht)^2;

if ind
%% Display en forma exponencial compleja
disp(['Hr = ', num2str(abs(Hr)), ' * exp(i*', num2str(angle(Hr)), ')']);
disp(['Ht = ', num2str(abs(Ht)), ' * exp(i*', num2str(angle(Ht)), ')']);
end

end


%[appendix]{"version":"1.0"}
%---
