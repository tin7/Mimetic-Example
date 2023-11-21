clear all;
addpath('C:\Users\x\Desktop\Metodos mimeticos\Example Mole\mole_MATLAB')
%%
% Valores de m y dt
valores_m = [10, 30, 60, 100, 500, 1000];


% Inicializar tabla de errores
errores = zeros(length(valores_m), 3); % Columnas para m, error Mimetico, error FR

for i = 1:length(valores_m)
    m = valores_m(i);
    % Spatial discretization
k = 2;         % Order of accuracy (spatial)
a = 0;         % Left boundary
b = 1;         % Right boundary
dx = (b-a)/m;  % Step length

verlet = 1;    % If verlet = 0 then it uses 4th order accurate for time (FR)

% 1D Staggered grid
xgrid = [a a+dx/2 : dx : b-dx/2 b];

% Mimetic operator (Laplacian)
L = lap(k, m, dx);
L(1, :) = 0;
L(end, :) = 0;
L(1, 1) = 1;
L(end, end) = 1;

% Wave propagation speed
c = 2;  % (T/p) Tension over density
% "Force" function
F = @(x) (c^2)*L*x;  % c^2 DivGrad x
% Simulation time
TIME = 0.5;

% Temporal discretization based on CFL condition
dt = dx/(4*c); % dt = h on Young's paper

% Initial condition
ICU = @(x) sin(pi*x) + sin(2*pi*x);      
ICV = @(x) zeros(m+1, 1); 


uold = ICU(xgrid');
vold = ICV(xgrid);
vold = [vold; vold(end)];

theta = 1/(2-2^(1/3)); % From Peter Young's paper

Nx = length(xgrid);
Nt = length(dt);
u_dif = zeros(Nx, Nt);
v_dif = zeros(size(u_dif));
CI= @(x) sin(pi*x) + sin(2*pi*x);
u_dif(:,1) = CI(xgrid');
u_dif(1, :) = 0; % u = 0 en x = a
u_dif(end, :) = 0; % u = 0 en x = b

u_exact = zeros(size(u_dif));

% Time integration loop
% v = VideoWriter('wave_simulation_1000.avi');
% open(v);


[X, T] = meshgrid(xgrid, 0:dt:TIME); % Crea matrices de malla para x y tiempo

unew_matrix = zeros(size(X));
udif_matrix = zeros(size(X));
uexact_matrix = zeros(size(X));

for t = 1 : TIME/dt
    % Apply "Verlet" algorithm (2nd-order in time)----------------
    if verlet
        uold = uold + 0.5*dt*vold;
        vnew = vold + dt*F(uold);
        unew = uold + 0.5*dt*vnew;
    % Apply "Forest-Ruth" algorithm (4th-order in time)--------------------
    else
        unew = uold + theta*0.5*dt*vold;
        vnew = vold + theta*dt*F(unew);
        unew = unew + (1-theta)*0.5*dt*vnew;
        vnew = vnew + (1-2*theta)*dt*F(unew);
        unew = unew + (1-theta)*0.5*dt*vnew;
        vnew = vnew + theta*dt*F(unew);
        unew = unew + theta*0.5*dt*vnew;
    end
    
    uold = unew;
    vold = vnew;
    %  Actualizar la posición usando la velocidad 
    %  de la mitad del paso de tiempo anterior  
    u_dif(:,t+1) = u_dif(:,t) + 0.5*dt*v_dif(:,t);
    %  Actualizar la velocidad usando la nueva posición
    for i = 2:Nx-1
        v_dif(i,t+1) = v_dif(i,t) + ...
        dt*(c^2/dx^2) * (u_dif(i+1,t+1) - 2*u_dif(i,t+1) + u_dif(i-1,t+1));
    end
    %  Actualizar la posición nuevamente usando la velocidad de la mitad
    %  del paso de tiempo siguiente
    u_dif(:,t+1) = u_dif(:,t+1) + 0.5*dt*v_dif(:,t+1);
    
    % Guardar los valores de los bordes
    u0_t(t+1) = u_dif(1, t+1);   % u(0, t)
    uL_t(t+1) = u_dif(end, t+1); % u(L, t)
    
    
    u_exact(:,t+1) = sin(pi*xgrid).*cos(2*pi*t*dt) + sin(2*pi*xgrid).*cos(4*pi*t*dt);
    drawnow;
    
    % Plot results
%     plot(xgrid, unew, 'LineWidth', 2)
%     title(['1D Wave equation \newlinet = ' num2str(dt*t)])
%     xlabel('x')
%     ylabel('u(x)')
%     axis([a b -1.5 1.5])
%     xticks([1, 2.5, 4])
%     xticklabels({'0', '0.5', '1'})
%     set(gcf, 'color', 'w')
%     grid on


% % Plot results
% plot(xgrid, unew, '-b', 'LineWidth', 2); % 'unew' en línea azul
% hold on; % Permite superponer las siguientes gráficas
% plot(xgrid, u_dif(:, t+1), '-r', 'LineWidth', 2); % 'u_dif' en línea roja
% hold on; % Permite superponer las siguientes gráficas
% plot(xgrid, u_exact(:,t+1), '-y', 'LineWidth', 2); 
% 
% title(['1D Wave equation \newlinet = ' num2str(dt*t)]);
% xlabel('x');
% ylabel('u(x)');
% axis([a b -1.5 1.5]);
% xticks([1, 2.5, 4]);
% xticklabels({'0', '0.5', '1'});
% set(gcf, 'color', 'w');
% grid on;
% legend('u\_new', 'u\_dif' , 'u\_exact'); % Agregar leyenda para distinguir las líneas
% hold off; % Desactivar superposición para gráficas futuras

unew_matrix(t+1,:) = unew;
udif_matrix(t+1,:) = u_dif(:,t+1);
uexact_matrix(t+1,:) = u_exact(:,t+1);


error_verlet = norm(unew_matrix - uexact_matrix, 'fro') / norm(uexact_matrix, 'fro');
error_fr = norm(udif_matrix - uexact_matrix, 'fro') / norm(uexact_matrix, 'fro');


    % Almacenar errores en la tabla
    errores(i, :) = [m, error_verlet, error_fr];


end


end

% Convertir a tabla para mejor visualización
tabla_errores = array2table(errores, 'VariableNames', {'m', 'Error_Verlet', 'Error_FR'});
% Filtrar las filas que no son todas ceros
tabla_filtrada = tabla_errores(any(tabla_errores{:,:} ~= 0, 2), :);

% Guardar la tabla filtrada en un archivo CSV
filename = 'errores_filtrados.csv';
writetable(tabla_filtrada, filename);
%%
% Graficar errores
figure; % Crea una nueva figura
loglog(errores(:,1), errores(:,2), '-o', 'LineWidth', 2); % Error Verlet
hold on; % Mantiene la figura actual para la siguiente grafica
loglog(errores(:,1), errores(:,3), '-x', 'LineWidth', 2); % Error FR

% Añadir detalles al gráfico
xlabel('m');
ylabel('Error');
title('Tasa de Convergencia vs m');
legend('Error Verlet', 'Error FD');
grid on;

% Guardar gráfico (opcional)
saveas(gcf, 'convergence_rate.png');

%%
%%% Graficar errores con plot y escala logarítmica
figure; % Crea una nueva figura
plot(errores(:,1), errores(:,2), '-o', 'LineWidth', 2); % Error Verlet
hold on; % Mantiene la figura actual para la siguiente grafica
plot(errores(:,1), errores(:,3), '-x', 'LineWidth', 2); % Error FR

% Establecer los ejes en escala logarítmica
set(gca, 'XScale', 'log', 'YScale', 'log');

% Añadir detalles al gráfico
xlabel('m');
ylabel('Error');
title('Tasa de Convergencia vs m');
legend('Error Mimetic', 'Error FD');
grid on;

% Guardar gráfico (opcional)
saveas(gcf, 'convergence_rate.png');

