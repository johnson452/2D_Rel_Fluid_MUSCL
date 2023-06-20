function [rho,ux,uy,uz,grid] = make_grid()

%%% Initialize memory %%%
%[DEFAULT] Setup Grid: (Boundary Grid):
grid.Nx = 51; % Only specified here
grid.Ny = 51;
grid.xmin = 0;
grid.xmax = 1.0;
grid.ymin = 0;
grid.ymax = 1.0;
grid.time = 0;
grid.t_max = 10.0;
grid.Output_interval = 20;
Nx = grid.Nx;
Ny = grid.Ny;

%[DEFAULT] Constants, updated in IC.m
grid.iter = 1;

%[DEFAULT] Grids, updated in IC.m
grid.dx = (grid.xmax - grid.xmin)/grid.Nx;
grid.dy = (grid.ymax - grid.ymin)/grid.Ny;
grid.time = 0;
grid.cfl = 0.98; %clf = udt/dx <= C_max
grid.dt = 0.98*sqrt(grid.dx*grid.dx + grid.dy*grid.dy)/10;
grid.NT = ceil(grid.t_max/grid.dt);
grid.Lx = (grid.xmax - grid.xmin);
grid.Ly = (grid.ymax - grid.ymin);

%Grid
grid.x = linspace(grid.xmin,grid.xmax,Nx);
grid.y = linspace(grid.ymin,grid.ymax,Ny);

%Build data
rho = zeros(Nx,Ny);
ux = zeros(Nx,Ny);
uy = zeros(Nx,Ny);
uz = zeros(Nx,Ny);


%Quantities
for j = 1:grid.Nx
    for k = 1:grid.Ny
        rho(j,k) = 2.5 + sin((2*grid.x(j)*pi)*((Nx-1)/Nx)/(grid.Lx))...
        + sin((2*grid.y(k)*pi)*((Ny-1)/Ny)/(grid.Ly));
        ux(j,k)  = 2.5 + sin((2*grid.x(j)*pi)*((Nx-1)/Nx)/(grid.Lx))...
        + sin((2*grid.y(k)*pi)*((Ny-1)/Ny)/(grid.Ly));
        uy(j,k)  = 2.5 + sin((2*grid.x(j)*pi)*((Nx-1)/Nx)/(grid.Lx))...
        + sin((2*grid.y(k)*pi)*((Ny-1)/Ny)/(grid.Ly));
        uz(j,k)  = 2.5 + sin((2*grid.x(j)*pi)*((Nx-1)/Nx)/(grid.Lx))...
        + sin((2*grid.y(k)*pi)*((Ny-1)/Ny)/(grid.Ly));
    end
end


%Right and left
grid.R = mod( linspace(1,Nx,Nx), Nx) + 1; %Good
grid.L = mod( linspace(-1,Nx-2,Nx), Nx) + 1; %Good

% KE vs t
grid.E_vs_t = zeros(1,grid.NT);
grid.time_vec = linspace(0,grid.t_max,grid.NT);

gamma = sqrt(1+ux.^2+uy.^2+uz.^2);
KE = (gamma - 1).*rho;
grid.E0 = sum(sum(KE)).*grid.dx;

%If the speed of light is violated by the IC
V2 = sqrt( (ux.^2+uy.^2+uz.^2)./(gamma.*gamma) );
if max(max(V2)) > 1
    fprintf("LIGHT SPEED VIOLATION in IC\n");
    exit();
end

end