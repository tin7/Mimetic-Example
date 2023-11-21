% Solves the 1D Wave equation with MDM and position Verlet algorithm
clc
clear all;
%%
addpath('../mole_MATLAB')

% Spatial discretization
k = 2;         % Order of accuracy (spatial)
m = 11;        % Number of cells
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


% Plot results
plot(xgrid, unew, '-b', 'LineWidth', 2); % 'unew' en línea azul
hold on; % Permite superponer las siguientes gráficas
plot(xgrid, u_dif(:, t+1), '-r', 'LineWidth', 2); % 'u_dif' en línea roja
hold on; % Permite superponer las siguientes gráficas
plot(xgrid, u_exact(:,t+1), '-y', 'LineWidth', 2); 

title(['1D Wave equation \newlinet = ' num2str(dt*t)]);
xlabel('x');
ylabel('u(x)');
axis([a b -1.5 1.5]);
xticks([1, 2.5, 4]);
xticklabels({'0', '0.5', '1'});
set(gcf, 'color', 'w');
grid on;
legend('u\_new', 'u\_dif' , 'u\_exact'); % Agregar leyenda para distinguir las líneas
hold off; % Desactivar superposición para gráficas futuras

unew_matrix(t+1,:) = unew;
udif_matrix(t+1,:) = u_dif(:,t+1);
uexact_matrix(t+1,:) = u_exact(:,t+1);


error_verlet = norm(unew_matrix - uexact_matrix, 'fro') / norm(uexact_matrix, 'fro');
error_fr = norm(udif_matrix - uexact_matrix, 'fro') / norm(uexact_matrix, 'fro');


end

% close(v)


%%
% Gráfica para unew_matrix
figure;
surf(X, T, unew_matrix);
title('unew');
xlabel('x');
ylabel('Time');
zlabel('u(x,t)');
shading interp;
colormap(jet);  
colorbar; 

% Gráfica para udif_matrix
figure;
surf(X, T, udif_matrix);
title('u\_dif');
xlabel('x');
ylabel('Time');
zlabel('u(x,t)');
shading interp;
colormap(jet);  
colorbar;  

% Gráfica para uexact_matrix
figure;
surf(X, T, uexact_matrix);
title('u\_exact');
xlabel('x');
ylabel('Time');
zlabel('u(x,t)');
shading interp;
colormap(jet);  
colorbar; 


%%
figure;
mesh(X, T, unew_matrix);  % Usar mesh en lugar de surf
title('unew');
xlabel('x');
ylabel('Time');
zlabel('u(x,t)');
colormap(jet);  % Mapa de colores jet
colorbar;  % Agrega una barra de colores para referencia

% Gráfica para udif_matrix
figure;
mesh(X, T, udif_matrix);  % Usar mesh en lugar de surf
title('u\_dif');
xlabel('x');
ylabel('Time');
zlabel('u(x,t)');
colormap(jet);  % Mapa de colores autumn
colorbar;  % Agrega una barra de colores para referencia

% Gráfica para uexact_matrix
figure;
mesh(X, T, uexact_matrix);  % Usar mesh en lugar de surf
title('u\_exact');
xlabel('x');
ylabel('Time');
zlabel('u(x,t)');
colormap(jet);  % Mapa de colores winter
colorbar;  % Agrega una barra de colores para referencia


%%
% Errores para unew respecto a uexact

% Error cuadrático medio (MSE)
MSE_unew = mean(mean((unew_matrix - uexact_matrix).^2));

% Error máximo relativo
MaxErrorRel_unew = max(max(abs((unew_matrix - uexact_matrix))));

% Norma L2
NormL2_unew = sqrt(sum(sum((unew_matrix - uexact_matrix).^2))*dx*dt);

% Errores para u_dif respecto a uexact

% Error cuadrático medio (MSE)
MSE_udif = mean(mean((udif_matrix - uexact_matrix).^2));

% Error máximo relativo
MaxErrorRel_udif = max(max(abs((udif_matrix - uexact_matrix))));

% Norma L2
NormL2_udif = sqrt(sum(sum((udif_matrix - uexact_matrix).^2))*dx*dt);

% Mostrando los resultados en la consola
fprintf('Para unew respecto a uexact:\n');
fprintf('MSE = %f\n', MSE_unew);
fprintf('Error Máximo  = %f\n', MaxErrorRel_unew);
fprintf('Norma L2 = %f\n\n', NormL2_unew);

fprintf('Para u_dif respecto a uexact:\n');
fprintf('MSE = %f\n', MSE_udif);
fprintf('Error Máximo  = %f\n', MaxErrorRel_udif);
fprintf('Norma L2 = %f\n', NormL2_udif);





