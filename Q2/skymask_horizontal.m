function skymask()
    % skymask - Creates a skymask visualization from a specific CSV file.
    
    % Load data from CSV
    data = readtable('skymask_A1_urban.csv');  % Specify your CSV file here
    
    % Extract azimuth and elevation angles
    azimuth = data.Azimuth_angle_deg;
    elevation = data.Elevation_angle_deg;

    % Create a figure with two subplots
    figure;

    %% Polar Plot
    subplot(1, 2, 1); % 1 row, 2 columns, first subplot
    azimuth_rad = deg2rad(azimuth); % Convert azimuth degrees to radians

    % Invert elevation for polar plot - center is 90, outer is 0
    polarplot(azimuth_rad, 90 - elevation, 'LineWidth', 2, 'Color', 'b');
    hold on;

    % Plot building outlines (example data, adjust as necessary)
    building_azimuth = [30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]; % Example data

    
    % Customize polar plot
    title('Skymask Visualization (Polar)');
    rlim([0 90]); % Limit radius to 90 degrees
    thetalim([0 360]); % Cover full azimuth range
    grid on;

    % Add legend
    legend('show', 'Location', 'bestoutside');
    
    %% Horizontal Plot
    subplot(1, 2, 2); % 1 row, 2 columns, second subplot
    plot(azimuth, elevation, 'LineWidth', 2, 'Color', 'b');
    title('Skymask Visualization (Horizontal)');
    xlabel('Azimuth Angle (degrees)');
    ylabel('Elevation Angle (degrees)');
    ylim([0 90]); % Limit elevation to 90 degrees
    grid on; % Display grid

    % Adjust layout
    sgtitle('Skymask Visualizations'); % Overall title for both subplots
end