%% Assignment 3 Part 2 - Andrew Paul
% This section of the assignmnet is taken from assignment 2 and plots the
% potnetial and electric field of the bottle-neck using the finite
% differences method

clear

nx = 50;
ny = 50;

% Create sparse G matrix
G = sparse(nx*ny,nx*ny);

% Conductivity outside box
sigma1 = 1;
% Conductivity inside box
sigma2 = 10^-2;

% Generate F matrix to set boundary conditions
F = zeros(nx*ny,1);

% Change for difference in bottle neck width
Lb = 0.4;
Wb = 0.6;

% Create matrix for mapping the conductivity and loop through to assign
% conductivity values for the given conditions
condMap = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        if (i>=Lb*nx && i<=Wb*nx && j<=Lb*ny) || (i>=Lb*nx && i<=Wb*nx && j>=Wb*ny)
            condMap(i,j) = sigma2;
        else
            condMap(i,j) = sigma1;
        end
    end
end

% Loop through to set boundary conditions

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 1;
            
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
            
        elseif j == 1
            nxm = j+(i-2)*ny;
            nxp = j+(i)*ny;
            nyp = j+1+(i-1)*ny;
            
            rxm = (condMap(i,j) + condMap(i-1,j))/2;
            rxp = (condMap(i,j) + condMap(i+1,j))/2;
            ryp = (condMap(i,j) + condMap(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
            
        elseif j == ny
            nxm = j+(i-2)*ny;
            nxp = j+(i)*ny;
            nym = j-1+(i-1)*ny;
            
            rxm = (condMap(i,j) + condMap(i-1,j))/2;
            rxp = (condMap(i,j) + condMap(i+1,j))/2;
            rym = (condMap(i,j) + condMap(i,j-1))/2;
        
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            
        else 
            nxm = j+(i-2)*ny;
            nxp = j+(i)*ny;
            nym = j-1+(i-1)*ny;
            nyp = j+1+(i-1)*ny;
            
            rxm = (condMap(i,j) + condMap(i-1,j))/2;
            rxp = (condMap(i,j) + condMap(i+1,j))/2;
            rym = (condMap(i,j) + condMap(i,j-1))/2;
            ryp = (condMap(i,j) + condMap(i,j+1))/2;
        
            G(n,n) = -(rxm+rxp+ryp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            
        end
        
    end
    
end

% Find voltage values using matrix operations
V = G\F;

% Create matrix to map voltage and loop through matrix to assign values
% from calculated voltage matrix
VMap = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        
        VMap(i,j) = V(n);
    end
end

% Plots

% Voltage plot
figure(2)
surf(VMap)
title('Voltage distribution over region')
xlabel('y')
ylabel('x')
colorbar
view(135,45)

%% 
% The voltage plot shown above gives the expected distribution as one side
% of the region has a voltage of V0 where the other side is set to zero.
% The barrier in the middle blocks uniform distribution to the other side
% of the region which is expected. As the conductivity inside the boxes is
% less than the conductivity outside the barrier boxes.

% Gradient used to plot electric field lines
[Ex,Ey] = gradient(VMap);

% Electric field plot
figure(3)
quiver(Ex,Ey)
title ('Electric Field Lines')
xlabel('x')
ylabel('y')
xlim([0 50])
ylim([0 50])

%% 
% The plot above shows that the electric field lines are strongest between
% the areas which have a lower conductivity. This is execpted as it is
% similar to the model of a parallel plate capcitor which has two plates of
% a larger conductivity seperated by a region with a lower conductivity
% creating a stronger electric field between the two higher conductivity
% regions.

