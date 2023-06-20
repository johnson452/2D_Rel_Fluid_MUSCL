function [rho,ux,uy,uz,grid] = push_once(rho,ux,uy,uz,grid)

% Construct Q and save for leapfrog
Q = construct(rho,ux,uy,uz, grid);
grid.Q_old = Q;

% Build quantitites
rho_u = zeros(3,grid.Nx,grid.Ny);
rho_u(1,:,:) = rho.*ux;
rho_u(2,:,:) = rho.*uy;
rho_u(3,:,:) = rho.*uz;

%Update quantities:
%Update old half quantities (n-1/2 -> n+1/2), Eqs 5, 6:
%[rho, rho_u] = y_push(rho, rho_u,...
%    rho, rho_u, grid);

%Set IC
rho0 = rho;
rho0_u0 = rho_u;

%SSP-RK3 (3 stage): push y
%Stage 1:
[rho, rho_u] = stage1y(rho, rho_u, grid);
%Stage 2:
[rho, rho_u] = stage2y(rho, rho_u, rho0, rho0_u0, grid);
%Stage 3:
[rho, rho_u] = stage3y(rho, rho_u, rho0, rho0_u0, grid);

%Set IC
rho0 = rho;
rho0_u0 = rho_u;

%SSP-RK3 (3 stage): push x
%Stage 1:
[rho, rho_u] = stage1x(rho, rho_u, grid);
%Stage 2:
[rho, rho_u] = stage2x(rho, rho_u, rho0, rho0_u0, grid);
%Stage 3:
[rho, rho_u] = stage3x(rho, rho_u, rho0, rho0_u0, grid);

% Update new quantities (n -> n + 1), Eqs 3, 4:
%[rho, rho_u] = x_push(rho, rho_u,...
%    rho, rho_u, grid);


%Retrieve u:
ux = squeeze(rho_u(1,:,:))./rho;
uy = squeeze(rho_u(2,:,:))./rho;
uz = squeeze(rho_u(3,:,:))./rho;

end


%construct Q
function Q = construct(N, Ux, Uy, Uz, grid)

% Address Q:
Q = zeros(4,grid.Nx,grid.Ny);

%Build Q : Q(i,:,:) = [ N ; N.*Ux ; N.*Uy ; N.*Uz ];
Q(1,:,:) = N;
Q(2,:,:) = N.*Ux;
Q(3,:,:) = N.*Uy;
Q(4,:,:) = N.*Uz;

end

%SSP-RK3 (stage 3)
function [rho, rho_u] = stage3x(rho, rho_u, rho0, rho0_u0, grid)

%Third stage calculation
[rho_star, rho_u_star] = x_push(rho, rho_u, grid);
rho = (1/3)*rho0 + (2/3)*rho_star;
rho_u = (1/3)*rho0_u0 + (2/3)*rho_u_star;

end

%SSP-RK3 (stage 2)
function [rho, rho_u] = stage2x(rho, rho_u, rho0, rho0_u0, grid)

%Second stage calculation
[rho_star, rho_u_star] = x_push(rho, rho_u, grid);
rho = (3/4)*rho0 + (1/4)*rho_star;
rho_u = (3/4)*rho0_u0 + (1/4)*rho_u_star;

end

%SSP-RK3 (stage 1)
function [rho, rho_u] = stage1x(rho, rho_u, grid)

%First stage calculation
[rho_star, rho_u_star] = x_push(rho, rho_u, grid);
rho = rho_star;
rho_u = rho_u_star;

end



%SSP-RK3 (stage 3)
function [rho, rho_u] = stage3y(rho, rho_u, rho0, rho0_u0, grid)

%Third stage calculation
[rho_star, rho_u_star] = y_push(rho, rho_u, grid);
rho = (1/3)*rho0 + (2/3)*rho_star;
rho_u = (1/3)*rho0_u0 + (2/3)*rho_u_star;

end

%SSP-RK3 (stage 2)
function [rho, rho_u] = stage2y(rho, rho_u, rho0, rho0_u0, grid)

%Second stage calculation
[rho_star, rho_u_star] = y_push(rho, rho_u, grid);
rho = (3/4)*rho0 + (1/4)*rho_star;
rho_u = (3/4)*rho0_u0 + (1/4)*rho_u_star;

end

%SSP-RK3 (stage 1)
function [rho, rho_u] = stage1y(rho, rho_u, grid)

%First stage calculation
[rho_star, rho_u_star] = y_push(rho, rho_u, grid);
rho = rho_star;
rho_u = rho_u_star;

end



%Leapfrog eqs 3, 4
function [rho, rho_u] = x_push(rho, rho_u, grid)


%Simplify the equations;
c = grid.dt/grid.dx;
R = grid.R;
L = grid.L;
ux = squeeze(rho_u(1,:,:))./rho;
uy = squeeze(rho_u(2,:,:))./rho;
uz = squeeze(rho_u(3,:,:))./rho;
gamma = sqrt(1 + ux.*ux + uy.*uy + uz.*uz);
vx = ux./gamma;
vy = uy./gamma;
vz = uz./gamma;

%Central averages to get ...
puvx_j_plus_half = zeros(3,grid.Nx,grid.Ny);
puvx_j_minus_half = zeros(3,grid.Nx,grid.Ny);
puvx_j_plus_half(1,:,:) = 0.5*(rho(R,:).*vx(R,:).*ux(R,:) + rho.*vx.*ux);
puvx_j_plus_half(2,:,:) = 0.5*(rho(R,:).*vx(R,:).*uy(R,:) + rho.*vx.*uy);
puvx_j_plus_half(3,:,:) = 0.5*(rho(R,:).*vx(R,:).*uz(R,:) + rho.*vx.*uz);
puvx_j_minus_half(1,:,:) = 0.5*(rho(L,:).*vx(L,:).*ux(L,:) + rho.*vx.*ux);
puvx_j_minus_half(2,:,:) = 0.5*(rho(L,:).*vx(L,:).*uy(L,:) + rho.*vx.*uy);
puvx_j_minus_half(3,:,:) = 0.5*(rho(L,:).*vx(L,:).*uz(L,:) + rho.*vx.*uz);

%Coef + or -
dv_j_plus = zeros(3,grid.Nx,grid.Ny);
dv_j_minus = zeros(3,grid.Nx,grid.Ny);
dv_j_plus(1,:,:) = (vx - vx(R,:));
dv_j_plus(2,:,:) = (vy - vy(R,:));
dv_j_plus(3,:,:) = (vz - vz(R,:));
dv_j_minus(1,:,:) = (vx(L,:) - vx);
dv_j_minus(2,:,:) = (vy(L,:) - vy);
dv_j_minus(3,:,:) = (vz(L,:) - vz);


%Coef + or -
coef_j_plus = zeros(3,grid.Nx,grid.Ny);
coef_j_minus = zeros(3,grid.Nx,grid.Ny);
for i = 1:3
coef_j_plus(i,:,:) = ( (gamma(R,:).*gamma)./(gamma(R,:)-gamma) ).*squeeze(dv_j_plus(i,:,:));
coef_j_minus(i,:,:) = ( (gamma.*gamma(L,:))./(gamma-gamma(L,:)) ).*squeeze(dv_j_minus(i,:,:));
end

%KEP constraint: (sum over all velocity directions)
pvx_j_plus_half = -squeeze(sum(coef_j_plus.*puvx_j_plus_half,1));
pvx_j_minus_half = -squeeze(sum(coef_j_minus.*puvx_j_minus_half,1));

%Update rho, u
rho = rho - c*(pvx_j_plus_half - pvx_j_minus_half);
rho_u = rho_u - c*(puvx_j_plus_half - puvx_j_minus_half);

end



%Leapfrog eqs 5, 6
function [rho_half, rho_u_half] = y_push(rho, rho_u, grid)


%Simplify the equations;
c = grid.dt/grid.dy;
R = grid.R;
L = grid.L;
ux = squeeze(rho_u(1,:,:))./rho;
uy = squeeze(rho_u(2,:,:))./rho;
uz = squeeze(rho_u(3,:,:))./rho;
gamma = sqrt(1 + ux.*ux + uy.*uy + uz.*uz);
vx = ux./gamma;
vy = uy./gamma;
vz = uz./gamma;

%Central averages to get ...
puvy_k_plus_half = zeros(3,grid.Nx,grid.Ny);
puvy_k_minus_half = zeros(3,grid.Nx,grid.Ny);
puvy_k_plus_half(1,:,:) = 0.5*(rho(:,R).*vy(:,R).*ux(:,R) + rho.*vy.*ux);
puvy_k_plus_half(2,:,:) = 0.5*(rho(:,R).*vy(:,R).*uy(:,R) + rho.*vy.*uy);
puvy_k_plus_half(3,:,:) = 0.5*(rho(:,R).*vy(:,R).*uz(:,R) + rho.*vy.*uz);
puvy_k_minus_half(1,:,:) = 0.5*(rho(:,L).*vy(:,L).*ux(:,L) + rho.*vy.*ux);
puvy_k_minus_half(2,:,:) = 0.5*(rho(:,L).*vy(:,L).*uy(:,L) + rho.*vy.*uy);
puvy_k_minus_half(3,:,:) = 0.5*(rho(:,L).*vy(:,L).*uz(:,L) + rho.*vy.*uz);

%Coef + or -
dv_k_plus = zeros(3,grid.Nx,grid.Ny);
dv_k_minus = zeros(3,grid.Nx,grid.Ny);
dv_k_plus(1,:,:) = (vx - vx(:,R));
dv_k_plus(2,:,:) = (vy - vy(:,R));
dv_k_plus(3,:,:) = (vz - vz(:,R));
dv_k_minus(1,:,:) = (vx(:,L) - vx);
dv_k_minus(2,:,:) = (vy(:,L) - vy);
dv_k_minus(3,:,:) = (vz(:,L) - vz);


%Coef + or -
coef_k_plus = zeros(3,grid.Nx,grid.Ny);
coef_k_minus = zeros(3,grid.Nx,grid.Ny);
for i = 1:3
coef_k_plus(i,:,:) = ( (gamma(:,R).*gamma)./(gamma(:,R)-gamma) ).*squeeze(dv_k_plus(i,:,:));
coef_k_minus(i,:,:) = ( (gamma.*gamma(:,L))./(gamma-gamma(:,L)) ).*squeeze(dv_k_minus(i,:,:));
end

%KEP constraint: (sum over all velocity directions)
pvy_k_plus_half = -squeeze(sum(coef_k_plus.*puvy_k_plus_half,1));
pvy_k_minus_half = -squeeze(sum(coef_k_minus.*puvy_k_minus_half,1));

%Update rho, u
rho_half = rho - c*(pvy_k_plus_half - pvy_k_minus_half);
rho_u_half = rho_u - c*(puvy_k_plus_half - puvy_k_minus_half);

end