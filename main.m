%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% 6/15/2023
% 2D MUSCL

%built-in-periodic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Main function
[rho,ux,uy,uz,grid] = make_grid();

% Push once
%[rho,ux,uy,uz,grid] = push_once(rho,ux,uy,uz,grid);

%Make the diagnostic Figure
figure('units','normalized','outerposition',[0 0 0.6 0.75])

%%% Time loop %%%
while(grid.time < grid.t_max)
    
    %Call i/o and diagnostics
    grid = diagnostics(rho,ux,uy,uz,grid);
    
    %Update the gridtime
    grid.time = grid.time + grid.dt;
    
    %Update the iterator
    grid.iter = grid.iter + 1;
 
    %Updater - updates all quantities simultaneosly
    % n -> n + 1 all quantities
    [rho,ux,uy,uz,grid] = push(rho,ux,uy,uz,grid);

end
%%% End Time Loop %%%
%%% End main %%%