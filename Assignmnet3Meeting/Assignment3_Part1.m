%% Assignment 3 Part 1 - Andrew Paul
% A constant voltage of 0.5 V was applied over the length of the
% semiconductor in the x direction (L) and the electric field (E) can be calculated
% by:
%
% $$E = V/L$$
%
% For a distance of 200 nm the electric field was calculated to be 0.5
% mV/nm or 500000 V/m.
%
% The force on each electron can be calculated by multiplying the charge of
% an electron (q) by the electric field in which the electron is present.
%
% $$F = Eq$$
%
% With an electron charge of 1.602 x $10^{-19}$ the force on the each
% electron was found to be 8.01 x $10^{-14}$ N.
%
% The acceleration of the electrons is calcaulted by dividing the force on
% the electrons by the effective mass of the electrons which is 0.26 x 9.11 x $10^{-31}$ kg.
%
% $$a = F/m$$
%
% The acceleration is calcualted to be 3.3817 x $10^{17}$.
%
% The formula for electron drift current density is given be the following
% formula:
%
% $$J_n = -e n V_d$$
%
% Where $J_n$ is the drift current density, $e$ is the charge of an electron,
% $n$ is the electron concentration (per $cm^2$)and $V_d$ is the average
% drift velocity.


%list of constants
m0 = 9.11e-31;
mn = 0.26*m0;
kB = 1.38e-23;
T = 300;
q = 1.602e-19;
n = 10^15; 

%region limits
xlim = 200e-9;
ylim = 100e-9;

%assignment 3 calcualtions
voltage = 0.1;
E = voltage/xlim;
F = E*q;
acc = F/mn;

%thermal velocity
vth = sqrt(2*kB*T/mn);

%initialize the number of electrons
num_electrons = 7;

%defining array for electrons (x postion, y position, angle, velocity)
electron = zeros(num_electrons, 6);

% the previous position of the electron (previous x position, previous y
% position)
electron_prev = zeros(num_electrons, 2);

%spacial step should be smaller than 1/100 of region size
time_step = xlim/vth/100;
time_total = time_step*250;
%num_step = time_total/time_step;

% used to make each electron a different colour
electron_colour = hsv(num_electrons);

% counter used to check temperature is constant
count = 0;

% scattering probability
Pscat = 1-exp(-time_step/0.2e-12);

%set an initial random postion and a fixed velocity for each electron
for i=1:num_electrons
    for j=1:6
        if(j==1)
            electron(i,j) = xlim*rand();
        elseif(j==2)
            electron(i,j) = ylim*rand();
        elseif(j==3)
            electron(i,j) = 2*pi*rand();
        elseif(j==4)
            electron(i,j) = randn()*vth;
        % set vx value
        elseif(j==5)
            electron(i,j) = electron(i,4)*cos(electron(i,3));
        % set vy value
        else
            electron(i,j) = electron(i,4)*sin(electron(i,3));
        end
    end
end

% figure(3)
% hist(electron(:,4))
% title('Velocity Distribution')

% define a temperature and time array for plotting
temperature= zeros(time_total/time_step,1);
time = zeros(time_total/time_step,1);

% counter for mean collision time
collision_count = 0;

running_time = 0;

% velocity array used to calculated mean free path
velocity = zeros(time_total/time_step,1);
avg_drift = zeros(time_total/time_step,1);
drift_current = zeros(time_total/time_step,1);
drift_vx = 0;

% update each electrons positon for each time step
for k=0:time_step:time_total
    avg_temp = 0;
    avg_velocity = 0;
    for m=1:num_electrons
        % allows electrons to pass through to the other side of the region
        %in the x-direction
        if (electron(m,1) >= xlim)
            electron(m,1) = 0;
            electron_prev(m,1) = 0;
        elseif (electron(m,1) <= 0)
            electron(m,1) = xlim;
            electron_prev(m,1) = xlim;
        end
        % electrons are reflected at the same angle if they strike the limits
        % of the region in the y-driection
        if ((electron(m,2) >= ylim) || (electron(m,2) <= 0))
            %electron(m,3) = pi - electron(m,3);
            %electron(m,4) = -electron(m,4);
            electron(m,6) = -electron(m,6);
        end
        
        % see if the particle scatters or not
        if(Pscat > rand())
            % scatters at a random angle
            electron(m,3) = 2*pi*rand();
            % new velocity for scattering - gaussian with some
            % MAXWELL-BOLTZMAN standard deviation
            vx_new = randn()*vth;
            vy_new = randn()*vth;
            v_new = sqrt(vx_new^2+vy_new^2);
            electron(m,4) = v_new;
            electron(m,5) = cos(electron(m,3))*v_new;
            electron(m,5) = sin(electron(m,3))*v_new;
            collision_count =+ 1;
        end
        
        avg_temp = avg_temp + (electron(m,4)^2)*mn/(2*kB);
        avg_velocity = avg_velocity + electron(m,4);
        
        %plot the movement of each electron
        if(k~=0)
            %UNCOMMENT TO SEE PLOT TRAJECTORY MOVIE
%             figure(1)
%             plot([electron_prev(m,1),electron(m,1)],[electron_prev(m,2),electron(m,2)],'color',electron_colour(m,:))
%             axis([0 xlim 0 ylim]);
        end
        
        drift_vx = drift_vx + electron(m,5);

    end
    
    %UNCOMMENT TO SEE PLOT TRAJECTORY MOVIE

%     title('Electron movement: random scattering')
%     xlabel('x-axis position (m)')
%     ylabel('y-axis position (m)')
%     hold on
%     pause(0.001)

   % set the previous postion of the electron to the current electron
   %postion for the next itteration
   electron_prev(:,1) = electron(:,1);
   electron_prev(:,2) = electron(:,2);
   
   % set the electron postion to an updated position
   electron(:,1) = electron(:,1) + electron(:,5).*time_step;
   electron(:,2) = electron(:,2) + electron(:,6).*time_step;
   
   electron(:,5) = electron(:,5) + acc*time_step;
   
   count = count +1;
   temperature(count,1) = avg_temp/num_electrons;
   time(count,1) = k + time_step;
   velocity(count,1) = avg_velocity;
   
   avg_drift(count,1) = drift_vx/num_electrons;
   
   drift_current(count,1) = q*n*avg_drift(count,1);

end

   figure(6)
   plot(time,drift_current)
   title('Drift current density over time')
   xlabel('time (s)')
   ylabel('Drift current density (A/cm^2)')
   
   figure(7)
   hist(drift_current')
   title('Drift Current Density Distribution')
   xlabel('Drift current density (A/cm^2)')
   ylabel('Counts')
   
% mean_collision = time_total/collision_count;
% avg_vth = 0;
% for n=1:500
%     avg_vth =+ velocity(n,1);
% end
% avg_vth = avg_vth/size(velocity,1);
% 
% MFP = avg_vth*mean_collision;
%    
% figure(2)
% plot(time,temperature)
% axis([0 time_total, 0 1100])
% title('Temperature of electrons over time')
% xlabel('time (s)')
% ylabel('Temperature (K)')

electron_grid = zeros(100,100);
temperature_grid = zeros(100,100);

% create density regions with grid vectors of final temperature and
%electron position
for x_pos=1:100
    for y_pos=1:100
        for q = 1:num_electrons
            if((electron(q,1) <= (xlim*(x_pos/100))) && (electron(q,1) > (xlim*((x_pos-1)/100))) && (electron(q,2) <= (ylim*(y_pos/100))) && (electron(q,2) > (ylim*((y_pos-1)/100))))
                electron_grid(x_pos,y_pos) =+ 1;
                temperature_grid(x_pos,y_pos) =+ (electron(q,4)^2)*mn/(2*kB);
            end
        end
    end
end


figure(2)
plot(time,temperature)
axis([0 time_total, 0 1100])
title('Temperature of electrons over time')
xlabel('time (s)')
ylabel('Temperature (K)')


figure(4)
surf(electron_grid)
colorbar
title('Electron Density Map')
xlabel('x-axis position (m)')
ylabel('y-axis position (m)')
shading interp
view(0,90)

% figure(5)
% surf(temperature_grid)
% colorbar
% title('Temperature Map')
% xlabel('x-axis position (m)')
% ylabel('y-axis position (m)')
% shading interp

%%
% 
% <<e_plot.png>>
% 


%% 
% The plot above shows the plot trajectory of 7 electrons. It is clear that
% there is a curve on the path of the electrons caused by an electric field
% applied across the region.
% 
% The plot trajectory shows that there is a curve on the path of the
% electrons due to the electric field applied across the region as
% expected. 
%
% Due to memory and timing concerns, the program was ran without the 
% plotting of the trajectories in order to proprly display the temperature 
% and dirft current density plots. 1000 electrons were used in the
% simulation and their results are shown below.
% 
%%
%
% <<d_curr_hist.PNG>>
% 
%%
%
% <<d_curr_plot.PNG>>
% 
%%
%
% <<e_density.PNG>>
%
%%
%
% <<temp_1.PNG>>
%
% 
% The plots for drift current density show that the initial current density
% is low which is expected because the initial velocity of the electrons is
% small and increases as the electic field accelerates the electrons. The
% histogram shows a relativly even distribution of the electrons which is
% due to the scattering which gives the electrons a new velocity, likely
% lower than the velocity they have been acceleratred to.
%
% The electron density map shows an relativly even distribution of
% electrons thoughout the area which is expected as there are no restricing
% spacial factors within the region.
%
% The temperature plot displays that the intial temperature of the
% electrons is lower and as they are accelerated and reach larger
% velocities, their temperature increases and somewhat saturates.



