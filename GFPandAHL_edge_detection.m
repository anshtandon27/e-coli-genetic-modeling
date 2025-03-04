INTEREST = 3;
%stran 1 =1, strain2 = 2. strain 3=3

% Plate radius (mm, call it Radius_Plate)
% Disk radius (mm, call it Radius_Disk)
% Experiment duration (minutes, call it T)
% Diffusion rate (mm^2/min, call it D)
% Source concentration (uM, call it sourceconc)

Radius_Plate = 85/2;
Radius_Disk = 3;
T = 21*60;

D = .15;
sourceconc = 10;
dAHL = 0;%4.8134e-4; %0;   %


%Optimized Model Params
param_names = {'LuxR', 'rho_R', 'delta_R', 'K_R', ...
               'alpha_TX_GFP', 'delta_TX_GFP', 'alpha_GFP', 'delta_GFP', 'n1'};

param_values = zeros(3, length(param_names));

% Strain 1 Parameters
strain_idx = 1;
param_values(strain_idx, :) = [2.47e-1, 0.5, 0.0231, 1.3e-3, 0.05, 0.2, 2.0, 4e-4, 1];

% Strain 2 Parameters
strain_idx = 2;
param_values(strain_idx, :) = [1e-2, 0.5, 0.0231, 1.3e-3, 0.05, 0.2, 2.0, 4e-4, 1];

% Strain 3 Parameters
strain_idx = 3;
param_values(strain_idx, :) = [1e-2, 0.5, 0.0231, 8.42e-5, 0.05, 0.2, 2.0, 4e-4, 1];

% Store model parameters separately
strainofinterest = param_values(INTEREST, :); %model_params; % param_values(#1-3, :); %%%%%%%%%%%%%%%%%%%%#####################################

% number of mesh points in the X- and Y-coordinates (use 201 points each to start)
% define the x- and y-step sizes (call them dx and dy). Remember MATLAB indices start at one
dx = 2* Radius_Plate/200;
dy = 2* Radius_Plate/200;

stability_factor = 0.2; % must be <=.25 for FTCS in 2D
dt = stability_factor*(dx^2)/D; % (min) time increment that fulfills the stability criterion (from Equation 6)
time = 0:dt:T; % Time Vector incremented in steps of dt

xgrid = -Radius_Plate:dx:Radius_Plate; % centers the grid in x around zero
ygrid = -Radius_Plate:dy:Radius_Plate; % centers the grid in y around zero
[X,Y] = meshgrid(xgrid, ygrid); % creates the matrix that describes the mesh using meshgrid

AHL_Initial = zeros(length(xgrid), length(ygrid));
% Set initial AHL concentration to zero across the whole plate described by the matrix AHL_Initial
[Disk_Indices_Row, Disk_Indices_Col] = find(sqrt(X.^2+Y.^2)<=Radius_Disk); 

% the loop below sets the initial concentration of the disk
for i = 1:length(Disk_Indices_Row)
    AHL_Initial(Disk_Indices_Row(i), Disk_Indices_Col(i)) = sourceconc;
end

% Define AHL matrix: (x points, y points, time points)
AHL = zeros(length(xgrid), length(ygrid), T);
% Initialize AHL at t=0 with AHL_Initial
AHL(:,:,1) = AHL_Initial;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate SUBSTRATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Fluxxing AHL
for t = 1:length(time) - 1  %time axis
    for i = 2:length(xgrid)-1  %x axis
        for j = 2:length(ygrid)-1 % y axis
            AHL(i,j,t+1) = AHL(i,j,t) +  dt * (D / dx^2) *(AHL(i-1, j, t) + AHL(i+1, j, t) + AHL(i, j-1, t) + AHL(i, j+1, t) - 4*AHL(i,j,t)) - dt*(dAHL*AHL(i,j,t));
        end
    end
    for i = 1:length(Disk_Indices_Row)
        AHL(Disk_Indices_Row(i), Disk_Indices_Col(i), t+1) = sourceconc;
    end
    % boundary conditions   
    AHL(1,:,t+1) = AHL(2,:,t); % Left boundary
    AHL(end,:,t+1) = AHL(end-1,:,t); % Right boundary
    AHL(:,1,t+1) = AHL(:,2,t); % Bottom boundary
    AHL(:,end,t+1) = AHL(:,end-1,t); % Top boundary
end

%Fluxxing R
R = zeros(length(xgrid), length(ygrid), length(time)); % Initial R = 0

LuxR = strainofinterest(1);  % LuxR concentration
rho_R = strainofinterest(2); % rho_R parameter
delta_R = strainofinterest(3); % delta_R parameter

for t = 1:length(time)-1  % Time axis
    for i = 2:length(xgrid)-1  % X axis
        for j = 2:length(ygrid)-1  % Y axis
            dR_dt = (rho_R * LuxR^2 * AHL(i,j,t)^2) - (delta_R * R(i,j,t));
            R(i,j,t+1) = R(i,j,t) + dt * dR_dt; 
        end
    end
end


%Flux TXGFP
K_R = strainofinterest(4); 
alpha_TXGFP = strainofinterest(5);
delta_TXGFP = strainofinterest(6);
n1 = strainofinterest(9);

TXGFP = zeros(size(R));

% Solve for TXGFP using the differential equation
for t = 1:length(time)-1  
    for i = 2:length(xgrid)-1 
        for j = 2:length(ygrid)-1 
            dTXGFP_dt = (alpha_TXGFP * (R(i,j,t)^n1) / (K_R^n1 + R(i,j,t)^n1)) - (delta_TXGFP * TXGFP(i,j,t));
            TXGFP(i,j,t+1) = TXGFP(i,j,t) + dt * dTXGFP_dt;
        end
    end
end

%Flux GFP
alpha_GFP = strainofinterest(7);
delta_GFP = strainofinterest(8);

GFP = zeros(size(TXGFP));

for t = 1:length(time)-1  
    for i = 2:length(xgrid)-1  
        for j = 2:length(ygrid)-1  
            dGFP_dt = (alpha_GFP * TXGFP(i,j,t)) - (delta_GFP * GFP(i,j,t));
            GFP(i,j,t+1) = GFP(i,j,t) + dt * dGFP_dt;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT SUBSTRATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%time points of interest in hours, .15 is used because when multipled by dt
%it = 1 which is the first point in the time index

time_indices = [dt, 3*60, 6*60, 9*60, 12*60, 15*60, 18*60, 21*60];

subplot_titles = ["t=0hrs", 't=3hrs', 't=6hrs', 't=9hrs', 't=12hrs','t=15hrs', 't=18hrs', 't=21hrs'];

% figure;
% % PLOT AHL
% for t = 1:length(time_indices)
%     subplot(2,4,t);
%     %remember to set the hour/min 0 timepoint to 1 index
%     allconc = AHL(:, :, round(time_indices(t)/dt));
%     imagesc(xgrid, ygrid, AHL(:,:,round(time_indices(t)/dt))); % Heat map
%     colorbar; % Add color scale
%     ylabel(colorbar, '[AHL] (uM)', 'Rotation', 270); % Set colorbar title parallel to the colorbar
%     title(subplot_titles(t));
%     axis equal;
%     xlabel('X (mm)');
%     ylabel('Y (mm)');
%     xlim([-Radius_Plate, Radius_Plate]);
%     ylim([-Radius_Plate, Radius_Plate]);
%     set(gcf, 'Position', [100, 100, 1200, 600]);
% end
% % 
% 
% %Plot R
% figure;
% for t =  1:length(time_indices)
%      subplot(2,4,t);
%     %remember to set the hour/min 0 timepoint to 1 index
%     if time_indices(t) == 1
%         time_indices(t) = 1;
%     end
%     allconc = R(:, :, round(time_indices(t)/dt));
%     imagesc(xgrid, ygrid, R(:,:,round(time_indices(t)/dt))); % Heat map
%     colorbar; % Add color scale
%     ylabel(colorbar, '[R] (uM)', 'Rotation', 270); % Set colorbar title parallel to the colorbar
%     title(subplot_titles(t));
%     axis equal;
%     xlabel('X (mm)');
%     ylabel('Y (mm)');
% 
%     xlim([-Radius_Plate, Radius_Plate]);
%     ylim([-Radius_Plate, Radius_Plate]);
%     set(gcf, 'Position', [100, 100, 1200, 600]);   
% end
% 
% 
% %Plot TXGFP
% figure;
% for t =  1:length(time_indices)
%      subplot(2,4,t);
%     %remember to set the hour/min 0 timepoint to 1 index
% 
%     allconc = TXGFP(:, :, round(time_indices(t)/dt));
%     imagesc(xgrid, ygrid, TXGFP(:,:,round(time_indices(t)/dt))); % Heat map
%     colorbar; % Add color scale
%     ylabel(colorbar, '[TXGFP] (uM)', 'Rotation', 270); % Set colorbar title parallel to the colorbar
%     title(subplot_titles(t));
%     axis equal;
%     xlabel('X (mm)');
%     ylabel('Y (mm)');
% 
%     xlim([-Radius_Plate, Radius_Plate]);
%     ylim([-Radius_Plate, Radius_Plate]);
%     set(gcf, 'Position', [100, 100, 1200, 600]);   
% end
% 
% %PLOT GFP
% 
% idx_21hr = find(time_indices == 21*60); % Assuming time is in minutes
% max_GFP_21hr = max(GFP(:, :, round(time_indices(idx_21hr)/dt)), [], 'all');
% % 
% figure;
% for t = 1:length(time_indices)
%     subplot(2, 4, t);
% 
%     % Extract GFP concentration at current time
%     allconc = GFP(:, :, round(time_indices(t)/dt));
% 
%     imagesc(xgrid, ygrid, allconc); % Heat map
%     colorbar; % Add color scale
%     ylabel(colorbar, '[GFP] (uM)', 'Rotation', 270);
%     title(subplot_titles(t));
%     axis equal;
%     xlabel('X (mm)');
%     ylabel('Y (mm)');
% 
%     xlim([-Radius_Plate, Radius_Plate]);
%     ylim([-Radius_Plate, Radius_Plate]);
% 
%     % Set the color scale limits for all plots
%     clim([0 max_GFP_21hr]); 
% 
%     set(gcf, 'Position', [100, 100, 1200, 600]);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%% GFP EDGE DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Threshold percentage (20% of the maximum GFP concentration)
%[.20, .0025, .0001]
threshold_percentage = .2;
radius_over_time = zeros(1, length(time)); 


for t = 1:length(time)
    % Extract valid indices outside the disk (y > Radius_Disk)
    valid_indices = find(ygrid > Radius_Disk); % Ignore values inside the disk

    % Extract centerline GFP values (assuming mid-column represents centerline)
    mid_x = floor(size(GFP, 2) / 2);  % Middle column index
    centerline_GFP(:, t) = GFP(valid_indices, mid_x, t); % Ensure it's a column vector

    % Find the maximum GFP value at the current time point (outside the disk)
    max_GFP = max(centerline_GFP(:, t));
    threshold_value = threshold_percentage * max_GFP;

    % Find indices where GFP concentration is above the threshold (excluding the disk)
    above_threshold_indices = find(centerline_GFP(:, t) >= threshold_value);
    
    if ~isempty(above_threshold_indices)
        % Define radius as the last index where GFP is above the threshold
        radius_index = above_threshold_indices(end);
        radius_over_time(t) = abs(ygrid(valid_indices(radius_index)));  % Update the radius
    else
        % If no values are above the threshold, set radius to the smallest valid y
        radius_over_time(t) = min(ygrid(valid_indices)); 
    end

    %edge case adjustments because model is physically restrained
    %center disk, edge border distance 
    if radius_over_time(t) < 3
       radius_over_time(t) = 0;
    end
    if radius_over_time(t) > 42.01
        radius_over_time(t) = 42.5;
    end
    if t < 180 && radius_over_time(t) > Radius_Plate / 2 - .01
        radius_over_time(t) = 0;  
    end
end


time_in_hours = time / (60);

figure;
plot(time_in_hours, radius_over_time, 'o-', 'LineWidth', 1, 'MarkerSize', 2, 'Color', "#2066a8");
hold on;

% Additional data (same as above)
VisDat1 = [1.013219228, 3.865196046, 11.24469162, 14.74483703, 17.50342293, 19.81348378, ...
           21.55539845, 23.39111547, 24.883198, 26.66797058, 27.96060285, 29.17842293, ...
           30.42045515, 31.18792537, 32.00624809, 32.75610571, 33.78482319, 34.68076498, ...
           35.64309462, 36.93721388, 38.27108061, 39.28051462];

VisE1 = [0.4775996483, 0.7037188979, 0.8825582517, 0.6724439488, 0.6792718855, 0.8496767292, ...
         0.9631014502, 0.8956070903, 0.9541937864, 1.050410215, 1.10886303, 1.111808414, ...
         1.242976952, 1.263666154, 1.284489237, 1.17762469, 0.9249270539, 0.8542195474, ...
         0.692516976, 0.60379789, 0.5521730118, 0.6467828667];


VisDat2 = [0.3575430516, 0.517395086, 1.7133236103, 2.489285602, 4.320985164, 6.438897179, ...
           10.18622731, 12.96294631, 14.59611128, 16.42780273, 17.86958304, 19.17336419, ...
           21.24826604, 22.97155682, 23.52735478, 24.83752838, 25.20520486, 27.96521035, ...
           27.12950597, 28.2140154, 29.01705743, 30.98915788];

VisE2 = [0.0229128341, 0.004117523577, 0.00808242226, 1.004115319, 1.403913899, 2.013771594, ...
         2.437105553, 1.710009991, 1.456419078, 0.9829097983, 0.8307950849, 0.5433046609, ...
         0.137376702, 0.002126053341, 0.06989728327, 0.01923097818, 0.4761124216, 0.6791991471, ...
         0.7131542409, 0.4595031126, 0.3411813757, 0.3581936621];


VisDat3 = [0.5463624946, 2.129158259, 3.626166503, 7.366374847, 9.810696408, 12.467208, ...
           14.38720758, 16.33766754, 18.38260237, 18.88561269, 20.70847908, 22.21292715, ...
           22.94696505, 24.08778358, 25.650265, 26.98960969, 28.0625153, 28.92563286, ...
           29.6054369, 29.75884059, 30.72239714, 30.77683822];

VisE3 = [0.2024216587, 0.3892275645, 0.02707489728, 0.03269421558, 0.04265573439, 0.05593775947, ...
         0.0651330076, 0.07100774947, 0.08122469183, 0.08684401014, 0.09348502268, 0.101403153, ...
         0.1052345064, 0.1121309425, 0.1213261906, 0.1279672032, 0.1346082157, 0.1389504162, ...
         0.142526346, 0.1445697345, 0.1471239701, 0.1494227821];


allData = zeros(3, length(VisDat3));
allData(1, :) = VisDat1;
allData(2, :) = VisDat2;
allData(3, :) = VisDat3;

allErrors = zeros(3, length(VisE3));
allErrors(1, :) = VisE1 *1.96;
allErrors(2, :) = VisE2*1.96;
allErrors(3, :) = VisE3*1.96;


dataofinterest = allData(INTEREST, :);
errorofinterest = allErrors(INTEREST, :);

% Time for the additional data in hours (same length as the additional data)
additional_time_in_hours = (0:length(dataofinterest)-1); 


% Initialize an array to store the predicted values from the original radius_over_time
predicted_values = NaN(size(dataofinterest));

% Find the nearest predicted values from radius_over_time based on additional_time_in_hours
for i = 1:length(additional_time_in_hours)
    % Find the closest time in the original time array (time_in_hours) to the current additional_time_in_hours(i)
    [~, idx] = min(abs(time_in_hours - additional_time_in_hours(i)));
    
    % Store the predicted value at that time
    predicted_values(i) = radius_over_time(idx);
end

squared_errors = (dataofinterest - predicted_values).^2;
range_radius = max(radius_over_time) - min(radius_over_time);
NMRSE = sqrt(mean(squared_errors)) / range_radius;




%%%%%%%%%%%%%%%%%%%%%%%%%%% AHL EDGE DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Threshold for normalized GFP (50% of maximum)
ahl_threshold = 0.5;

% Calculate the AHL threshold concentration for the strain of interest
AHL_conc = logspace(-4, 1, 100); % AHL concentration range from 0.0001 to 10 μM

% Extract parameters for strain of interest
strain_params = param_values(INTEREST, :);

% Extract individual parameters
LuxR = strain_params(1);
rho_R = strain_params(2);
delta_R = strain_params(3);
K_R = strain_params(4);
alpha_TXGFP = strain_params(5);
delta_TXGFP = strain_params(6);
alpha_GFP = strain_params(7);
delta_GFP = strain_params(8);
n1 = strain_params(9);

% Calculate R, TXGFP, and GFP for each AHL concentration
ahlR = (rho_R * LuxR^2 * AHL_conc.^2) / delta_R;
ahlTXGFP = (alpha_TXGFP * (ahlR.^n1) ./ (K_R^n1 + ahlR.^n1)) / delta_TXGFP;
ahlGFP = (alpha_GFP * ahlTXGFP) / delta_GFP;

% Normalize GFP
normalized_GFP = ahlGFP / max(ahlGFP);

% Find AHL concentration at 50% max GFP (threshold)
[~, idx] = min(abs(normalized_GFP - ahl_threshold));
ahl_threshold_conc = AHL_conc(idx);

fprintf('Strain %d: AHL threshold concentration = %.6f μM\n', INTEREST, ahl_threshold_conc);

% Initialize array to store AHL edge radius over time
ahl_radius_over_time = zeros(1, length(time));

% Track where the AHL threshold concentration occurs in the spatial model
for t = 1:length(time)
    % Extract valid indices outside the disk (y > Radius_Disk)
    valid_indices = find(ygrid > Radius_Disk); % Ignore values inside the disk

    % Extract centerline AHL values
    mid_x = floor(size(AHL, 2) / 2);  % Middle column index
    centerline_AHL = AHL(valid_indices, mid_x, t);
    
    % Find where AHL concentration is closest to the threshold
    [~, idx] = min(abs(centerline_AHL - ahl_threshold_conc));
    
    if ~isempty(idx)
        % Store the distance from the origin (radius)
        ahl_radius_over_time(t) = abs(ygrid(valid_indices(idx)));
    else
        ahl_radius_over_time(t) = 0;
    end


    if ahl_radius_over_time(t) < 3
        ahl_radius_over_time(t) = 0;
    end

    if t < 180 && ahl_radius_over_time(t) > Radius_Plate / 2 - 0.01 - Radius_Disk
        ahl_radius_over_time(t) = 0;
    end
end



% Plot AHL edge
plot(time_in_hours, ahl_radius_over_time, 's-', 'LineWidth', 1, 'MarkerSize', 2, 'Color', '#9f8c8c');

% Plot experimental data points
errorbar(additional_time_in_hours, dataofinterest, errorofinterest, 'o', 'LineStyle', 'none', 'LineWidth', 1, 'MarkerFaceColor', '#ae282c', 'MarkerEdgeColor', '#ae282c', 'Color', '#ae282c', 'MarkerSize', 3);


% Add labels, title, and grid
xlabel('Time (hours)');
ylabel('Radius (mm)');
ylim([0 45]);

title(sprintf('GFP and AHL Edge VS Time (S%d, %.0f%% Threshold)', INTEREST, threshold_percentage * 100));

% Legend text with NMRSE
legend_text = sprintf('NMRSE: %.4f', NMRSE);

% Add an invisible plot to allow NMRSE to show in the legend
plot(NaN, NaN, 's', 'Color', 'none'); % This adds an invisible entry with no icon

% Add the legend with Model Prediction, AHL Edge, Visual Data, and NMRSE
legend({'Model GFP Edge', 'Model AHL Edge', 'GFP Photo Edge', legend_text}, 'Location', 'southeast');

grid on;
hold off;

