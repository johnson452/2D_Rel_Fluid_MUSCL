%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/17/2023
%Fluid algroithm, V (divergence U)

%Notes:
%-2D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Update the quanitites Ux, Uy, Uz (t -> t + dt)
function [N,Ux,Uy,Uz,grid] = push(N,Ux,Uy,Uz,grid)

% Build Q
Q = construct(N, Ux, Uy, Uz, grid);
Nx = grid.Nx;
Ny = grid.Ny;

% Grid R and L:
R = grid.R;
L = grid.L;

%Reconstruct U(x) (soln) within one cell, also slope dU
dQx = reconstruct(Q,Q(:,L,:),Q(:,R,:));
dQy = reconstruct(Q,Q(:,:,L),Q(:,:,R));

%Update solution, primative variablesabs
%Q_tilde = Q - grid.dt/(2*grid.dx)*AQ(Q,grid).*dQ;
Ax = AQx(Q,grid);
Ay = AQy(Q,grid);
Q_tilde = zeros(4,Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        Q_tilde(:,i,j) = Q(:,i,j) - grid.dt/(2*grid.dx)*( Ax(:,:,i,j)*dQx(:,i,j) + Ay(:,:,i,j)*dQy(:,i,j));
    end
end


%Matrix comps
% for i = 1:4
%     for j = 1:4
%         subplot(2,1,1)
%         surf(diff(squeeze(Ax(i,j,:,:)),1))
%         subplot(2,1,2)
%         surf(diff(squeeze(Ay(i,j,:,:)),1))
%         fprintf("i = %d, j = %d\n",i,j)
%         pause(5)
%         clf()
%     end
% end


%Get edge values (edges_linear)
% Output: (cell edge values) in cell i
% [ i - 1/2 (+), i + 1/2 (-) ]
[Q_plus_I_x, Q_minus_I_x] = edges_linear(Q_tilde,dQx);
[Q_plus_R_x, ~] = edges_linear(Q_tilde(:,R,:),dQx(:,R,:));
[~, Q_minus_L_x] = edges_linear(Q_tilde(:,L,:),dQx(:,L,:));
[Q_plus_I_y, Q_minus_I_y] = edges_linear(Q_tilde,dQy);
[Q_plus_R_y, ~] = edges_linear(Q_tilde(:,:,R),dQy(:,:,R));
[~, Q_minus_L_y] = edges_linear(Q_tilde(:,:,L),dQy(:,:,L));

%Get the edge values (cell edges, in cell i)
% [Q_plus_I, Q_minus_I] = edges(Q_tilde);
% [Q_plus_R, ~] = edges(Q_tilde(:,R));
% [~, Q_minus_L] = edges(Q_tilde(:,L));

%Update the fluxes
F_R_x =  Flux_x(Q_minus_I_x,Q_plus_R_x,grid);
F_L_x = Flux_x(Q_minus_L_x,Q_plus_I_x,grid);
F_R_y =  Flux_y(Q_minus_I_y,Q_plus_R_y,grid);
F_L_y = Flux_y(Q_minus_L_y,Q_plus_I_y,grid);

%Compute the updated Q
Q = Q - grid.dt/(grid.dx)*(F_R_x - F_L_x) ...
    - grid.dt/(grid.dy)*(F_R_y - F_L_y);

% Destruct Q into it's components
[N, Ux, Uy, Uz] = destruct(Q);

end


%Locally defined functions for MUSCL-Handcock:


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


%Destruct Q into N, Ux, Uy, Uz
function [N, Ux, Uy, Uz] = destruct(Q)
N =  squeeze(Q(1,:,:));
Ux = squeeze(Q(2,:,:))./N;
Uy = squeeze(Q(3,:,:))./N;
Uz = squeeze(Q(4,:,:))./N;
end

%Fluxes N
function [Fl] = Flux_x(QL,QR,grid)

[~, UxR, UyR, UzR] = destruct(QR);
[~, UxL, UyL, UzL] = destruct(QL);

% compute c for the lax flux
UR2 = ( UxR.^2 + UyR.^2 + UzR.^2 );
UL2 = ( UxL.^2 + UyL.^2 + UzL.^2 );
vRx = UxR./(sqrt(1 + UR2));
vLx = UxL./(sqrt(1 + UL2));
c = max( abs(vRx), abs(vLx) ); %%% FIX ABS | lamdba^p|

%Rusanov Flux F = vxQ
F_cQ = zeros(4,grid.Nx,grid.Ny);
for i = 1:4
    F_cQ(i,:,:) = (1/2) *(  vRx.*squeeze(QR(i,:,:)) + vLx.*squeeze(QL(i,:,:)) ...
        - c.*squeeze( QR(i,:,:) - QL(i,:,:) ) );
end
Fl =  ( F_cQ  );

end

%Fluxes N
function [Fl] = Flux_y(QL,QR,grid)

[~, UxR, UyR, UzR] = destruct(QR);
[~, UxL, UyL, UzL] = destruct(QL);

% compute c for the lax flux
UR2 = ( UxR.^2 + UyR.^2 + UzR.^2 );
UL2 = ( UxL.^2 + UyL.^2 + UzL.^2 );
vRy = UyR./(sqrt(1 + UR2));
vLy = UyL./(sqrt(1 + UL2));
c = max( abs(vRy), abs(vLy) ); %%% FIX ABS | lamdba^p|

%Rusanov Flux F = vxQ
F_cQ = zeros(4,grid.Nx,grid.Ny);
for i = 1:4
    F_cQ(i,:,:) = (1/2) * ( vRy.*squeeze(QR(i,:,:)) + vLy.*squeeze(QL(i,:,:)) ...
        - c.*squeeze( QR(i,:,:) - QL(i,:,:) ) );
end
Fl =  ( F_cQ );

end

%Reconstruction
function [dW] = reconstruct(Wi,Wm,Wp)

%Average Dw
dW = ave( Wi - Wm, Wp - Wi );

end

% Averaging
function [dW] = ave( Wm, Wp )

avg_type = "minmod"; %"Supebee"; %"standard";
a = Wm; b = Wp;
ab = a.*b;
sz_a = size(a);
dW = zeros(sz_a);

% Standard Averaging
if avg_type == "standard"
    dW = (Wm + Wp)/2;
elseif avg_type == "minmod"
    for i = 1:sz_a(1)
        for j = 1:sz_a(2)
            for k = 1:sz_a(3)
                if ab(i,j,k) > 0
                    dW(i,j,k) = minmod( [(a(i,j,k) + b(i,j,k))/2 , 2*a(i,j,k), 2*b(i,j,k)] );
                else
                    dW(i,j,k) = 0;
                end
            end
        end
    end
elseif avg_type == "Supebee"
    for i = 1:sz_a(1)
        for j = 1:sz_a(2)
            for k = 1:sz_a(3)
                if ab(i,j,k) > 0
                    max_v =  maxmod(  [ a(i,j,k)  , b(i,j,k)  ] );
                    min_v = minmod (  [ 2*a(i,j,k), 2*b(i,j,k)] );
                    dW(i,j,k) = minmod( [ max_v   , min_v   ] );
                else
                    dW(i,j,k) = 0;
                end
            end
        end
    end
end


end

% Minmod used for averaging
function [val] = minmod(a)
if (max(a) > 0) && (min(a) >  0)
    val = min(a);
elseif (max(a) < 0) && (min(a) <  0)
    val = max(a);
else
    val = 0;
end
end

%Maxmod used for averaging
function [val] = maxmod(a)
if (max(a) > 0) && (min(a) >  0)
    val = max(a);
elseif (max(a) < 0) && (min(a) <  0)
    val = min(a);
else
    val = 0;
end
end


% Function Ax(Q)
function [A] = AQx(Q,grid)

% Build the 4x4 flux Jacobian... (See Maxima generated values)
% [ A11, A12, A13, A14 ]
% [ A21, A22, A23, A24 ]
% [ A31, A32, A33, A34 ]
% [ A41, A42, A43, A44 ]

%Recover primitive variables
[~, Ux, Uy, Uz] = destruct(Q);

%Definte helper functions
a = Uz.^4+(2*Uy.^2+2*Ux.^2+2).*Uz.^2+Uy.^4+(2*Ux.^2+2).*Uy.^2+Ux.^4+2*Ux.^2+1;
gamma = sqrt(1 + Ux.*Ux + Uy.*Uy + Uz.*Uz );

%Build the components of A
A11 = ((Ux.*Uz.^2+Ux.*Uy.^2+Ux.^3).*gamma)./a;
A12 = ((Uz.^2+Uy.^2+1).*gamma)./a;
A13 = -(Ux.*Uy)./gamma.^3;
A14 = -(Ux.*Uz)./gamma.^3;
A21 = -Ux.^2./gamma.^3;
A22 = ((2.*Ux.*Uz.^2+2.*Ux.*Uy.^2+Ux.^3+2.*Ux).*gamma)./a;
A23 = -(Ux.^2.*Uy)./gamma.^3;
A24 = -(Ux.^2.*Uz)./gamma.^3;
A31 = -(Ux.*Uy)./gamma.^3;
A32 = ((Uy.*Uz.^2+Uy.^3+Uy).*gamma)./a;
A33 = ((Ux.*Uz.^2+Ux.^3+Ux).*gamma)./a;
A34 = -(Ux.*Uy.*Uz)./gamma.^3;
A41 = -(Ux.*Uz)./gamma.^3;
A42 = ((Uz.^3+(Uy.^2+1).*Uz).*gamma)./a;
A43 = -(Ux.*Uy.*Uz)./gamma.^3;
A44 = ((Ux.*Uy.^2+Ux.^3+Ux).*gamma)./a;

%Assemble A
A = zeros(4,4,grid.Nx,grid.Ny);
for i = 1:grid.Nx
    for j = 1:grid.Ny
        A(:,:,i,j) = [ [ A11(i,j), A12(i,j), A13(i,j), A14(i,j) ];...
            [ A21(i,j), A22(i,j), A23(i,j), A24(i,j) ];...
            [ A31(i,j), A32(i,j), A33(i,j), A34(i,j) ];...
            [ A41(i,j), A42(i,j), A43(i,j), A44(i,j) ] ];
    end
end

end




% Function Ax(Q)
function [A] = AQy(Q,grid)

% Build the 4x4 flux Jacobian... (See Maxima generated values)
% [ A11, A12, A13, A14 ]
% [ A21, A22, A23, A24 ]
% [ A31, A32, A33, A34 ]
% [ A41, A42, A43, A44 ]

%Recover primitive variables
[~, Ux, Uy, Uz] = destruct(Q);

%Definte helper functions
a = Uz.^4+(2*Uy.^2+2*Ux.^2+2).*Uz.^2+Uy.^4+(2*Ux.^2+2).*Uy.^2+Ux.^4+2*Ux.^2+1;
gamma = sqrt(1 + Ux.*Ux + Uy.*Uy + Uz.*Uz );

%Build the components of A
A11 = ((Uy.*Uz.^2+Uy.^3+Ux.^2.*Uy).*gamma)./a;
A12 = -(Ux.*Uy)./gamma.^3;
A13 = ((Uz.^2+Ux.^2+1).*gamma)./a;
A14 = -(Uy.*Uz)./gamma.^3;

A21 = -(Ux.*Uy)./gamma.^3;
A22 = ((Uy.*Uz.^2+Uy.^3+Uy).*gamma)./a;
A23 = ((Ux.*Uz.^2+Ux.^3+Ux).*gamma)./a;
A24 = -(Ux.*Uy.*Uz)./gamma.^3;

A31 = -Uy.^2./gamma.^3;
A32 = -(Ux.*Uy.^2)./gamma.^3;
A33 = ((2*Uy.*Uz.^2+Uy.^3+(2*Ux.^2+2).*Uy).*gamma)./a;
A34 = -(Uy.^2.*Uz)./gamma.^3;

A41 = -(Uy.*Uz)./gamma.^3;
A42 = -(Ux.*Uy.*Uz)./gamma.^3;
A43 = ((Uz.^3+(Ux.^2+1).*Uz).*gamma)./a;
A44 = ((Uy.^3+(Ux.^2+1).*Uy).*gamma)./a;

%Assemble A
A = zeros(4,4,grid.Nx,grid.Ny);
for i = 1:grid.Nx
    for j = 1:grid.Ny
        A(:,:,i,j) = [ [ A11(i,j), A12(i,j), A13(i,j), A14(i,j) ];...
            [ A21(i,j), A22(i,j), A23(i,j), A24(i,j) ];...
            [ A31(i,j), A32(i,j), A33(i,j), A34(i,j) ];...
            [ A41(i,j), A42(i,j), A43(i,j), A44(i,j) ] ];
    end
end

end

function [W_plus, W_minus] = edges_linear(W_tilde,dW)
W_plus = W_tilde - dW/2;
W_minus = W_tilde + dW/2;
end
