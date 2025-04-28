function skymask_RAIM_enhanced()
    % Enhanced RAIM implementation with improved visualization
    % This function loads navigation solution data, performs RAIM calculations,
    % and creates visualizations of the selected results
    
    %% Load and verify data
    filePath = 'navSolution_OpenSky.mat';
    try
        navData = load(filePath);
        pseudoranges = navData.navSolutions.correctedP;         % 5×178 matrix
        satellite_positions = navData.navSolutions.satPositions; % 3×5×178 matrix
        
        % Verify data dimensions
        [numSatellites, numEpochs] = size(pseudoranges);
        fprintf('Data loaded: %d satellites across %d epochs\n', numSatellites, numEpochs);
        
        if size(satellite_positions, 1) ~= 3 || size(satellite_positions, 2) ~= numSatellites || size(satellite_positions, 3) ~= numEpochs
            error('Mismatched dimensions between pseudoranges and satellite positions');
        end
    catch ME
        error('Error loading navigation data: %s', ME.message);
    end
    
    %% Configuration parameters
    config.sigma = 3;                     % Standard deviation of pseudorange measurements (meters)
    config.alpha = 0.02;                  % Significance level for chi-square test
    config.p_md = 1e-7;                   % Probability of missed detection
    config.criticalValue = chi2inv(1-config.alpha, numSatellites-4);
    config.k = chi2inv(1-config.p_md, 1); % For Protection Level calculation
    config.alert_limit = 50;              % Horizontal Alert Limit (HAL) in meters
    
    %% Initialize result arrays
    results.positions = zeros(numEpochs, 4);       % [x, y, z, clock_bias]
    results.residuals = zeros(numSatellites, numEpochs);
    results.chi_square = zeros(numEpochs, 1);
    results.protectionLevels = zeros(numEpochs, 1);
    results.faultDetected = false(numEpochs, 1);
    
    % Define a "true" reference position for error calculation
    reference_position = [0, 0, 0];
    
    %% Process each epoch
    for epoch = 1:numEpochs
        % Extract data for current epoch
        current_pseudoranges = pseudoranges(:, epoch);
        
        % Reshape satellite positions to 5×3 matrix (satellites × coordinates)
        current_sat_positions = reshape(satellite_positions(:, :, epoch), 3, numSatellites)';
        
        % Skip epoch if data is invalid
        if any(isnan(current_sat_positions(:))) || any(isnan(current_pseudoranges))
            fprintf('Warning: Missing data at epoch %d - skipping\n', epoch);
            continue;
        end
        
        % Build design matrix A
        A = zeros(numSatellites, 4);
        for i = 1:numSatellites
            sat_pos = current_sat_positions(i, :);
            sat_range = norm(sat_pos);
            
            % Avoid division by zero
            if sat_range > 1e-10
                A(i, 1:3) = sat_pos / sat_range;
            else
                warning('Very small satellite range at epoch %d, satellite %d', epoch, i);
                A(i, 1:3) = [0, 0, 0];
            end
            
            A(i, 4) = -1; % Clock bias component
        end
        
        % Initial weights (identity matrix assuming equal weights)
        W = eye(numSatellites);
        
        % Compute position estimate using Weighted Least Squares
        try
            % Add small regularization for numerical stability
            position = (A' * W * A + 1e-10 * eye(4)) \ (A' * W * current_pseudoranges);
            
            % Compute measurement residuals
            residuals = current_pseudoranges - A * position;
            
            % Variance of residuals
            sigma_r2 = (residuals' * residuals) / (numSatellites - 4);
            
            % Chi-square test statistic
            chi_square = (residuals' * W * residuals) / sigma_r2;
            
            % Fault detection based on chi-square test
            fault_detected = chi_square > config.criticalValue;
            
            % Position error covariance matrix
            P = inv(A' * W * A);
            
            % Store results
            results.positions(epoch, :) = position';
            results.residuals(:, epoch) = residuals;
            results.chi_square(epoch) = chi_square;
            results.protectionLevels(epoch) = config.k * config.sigma; % Original PL
            results.faultDetected(epoch) = fault_detected;
            
            if fault_detected
                fprintf('Epoch %d: Fault detected (chi-square = %.2f > %.2f)\n', ...
                    epoch, chi_square, config.criticalValue);
            end
        catch ME
            warning('Computation error at epoch %d: %s', epoch, ME.message);
        end
    end
    
    %% Compute a reference position as the mean of all valid position estimates
    valid_positions = results.positions(sum(abs(results.positions), 2) > 0, 1:3);
    reference_position = mean(valid_positions, 1);
    
    %% Create visualizations - now with four plots
    visualizeResults(satellite_positions, results, numEpochs, config);
end

function visualizeResults(satellite_positions, results, numEpochs, config)
    % Create visualization of RAIM results with four plots
    
    % Create figure with subplots - 2x2 grid for four plots
    figure('Position', [100, 100, 1200, 900], 'Color', 'white');
    
    % Plot 1: 3D satellite constellation
    ax1 = subplot(2, 2, 1);
    plotSatelliteConstellation3D(satellite_positions, numEpochs, ax1);
    
    % Plot 2: XY view of satellite constellation (new)
    ax2 = subplot(2, 2, 2);
    plotSatelliteConstellationXY(satellite_positions, numEpochs, ax2);
    
    % Plot 3: Chi-square test statistics over time
    ax3 = subplot(2, 2, 3);
    plotChiSquareStatistics(results.chi_square, config.criticalValue, numEpochs, ax3);
    
    % Plot 4: Protection levels over time
    ax4 = subplot(2, 2, 4);
    plotProtectionLevels(results.protectionLevels, numEpochs, ax4);
    
    % Add overall title and adjust spacing
    sgtitle('RAIM Analysis Results', 'FontSize', 14, 'FontWeight', 'bold');
    set(gcf, 'Color', 'white');
    
    % Improve figure appearance
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 1.5);
    
    % Improve layout
    tight_subplot = false;
    if exist('tight_subplot', 'file') == 2 && tight_subplot
        tight_subplot(2, 2, [0.08 0.03], [0.08 0.05], [0.05 0.02]);
    end
end

function plotSatelliteConstellation3D(satellite_positions, numEpochs, ax)
    % Visualize satellite constellation in 3D with a modern style
    
    % Select subset of epochs to avoid overcrowding (every 10th epoch)
    epoch_subset = 1:10:numEpochs;
    
    % Define colormap
    cmap = turbo(length(epoch_subset));
    
    % Create axes
    axes(ax);
    hold on;
    
    % Plot Earth (approximate size)
    [X, Y, Z] = sphere(50);
    scale_factor = 1e7; % Scale
    earth_handle = surf(X*scale_factor, Y*scale_factor, Z*scale_factor, ...
        'FaceColor', [0.3, 0.5, 0.9], 'FaceAlpha', 0.3, ...
        'EdgeColor', [0.4, 0.6, 1], 'EdgeAlpha', 0.2);
    
    % Plot satellites for selected epochs
    for i = 1:length(epoch_subset)
        epoch = epoch_subset(i);
        current_sat_positions = reshape(satellite_positions(:, :, epoch), 3, [])';
        
        % Plot satellites
        scatter3(current_sat_positions(:, 1), current_sat_positions(:, 2), ...
            current_sat_positions(:, 3), 80, cmap(i,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end
    
    % Add colorbar
    c = colorbar;
    c.Label.String = 'Epoch Number';
    clim([min(epoch_subset), max(epoch_subset)]);
    colormap(ax, turbo);
    
    % Set axis properties
    grid on;
    box on;
    axis equal;
    
    % Format axes with appropriate scaling
    ax.XAxis.Exponent = 7;
    ax.YAxis.Exponent = 7;
    ax.ZAxis.Exponent = 7;
    
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('3D Satellite Constellation');
    
    % Set view angle
    view(45, 30);
    
    % Add lighting
    light('Position', [1 1 1], 'Style', 'infinite');
    lighting phong;
    
    hold off;
end

function plotSatelliteConstellationXY(satellite_positions, numEpochs, ax)
    % Visualize satellite constellation in 2D (X-Y view)
    
    % Select subset of epochs to avoid overcrowding (every 10th epoch)
    epoch_subset = 1:10:numEpochs;
    
    % Define colormap
    cmap = turbo(length(epoch_subset));
    
    % Create axes
    axes(ax);
    hold on;
    
    % Plot Earth (approximate size)
    theta = linspace(0, 2*pi, 100);
    scale_factor = 1e7; % Scale
    earth_x = scale_factor * cos(theta);
    earth_y = scale_factor * sin(theta);
    fill(earth_x, earth_y, [0.3, 0.5, 0.9], 'FaceAlpha', 0.3, 'EdgeColor', [0.4, 0.6, 1]);
    
    % Plot satellites for selected epochs
    for i = 1:length(epoch_subset)
        epoch = epoch_subset(i);
        current_sat_positions = reshape(satellite_positions(:, :, epoch), 3, [])';
        
        % Plot satellites (X-Y projection)
        scatter(current_sat_positions(:, 1), current_sat_positions(:, 2), 100, ...
            cmap(i,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end
    
    % Add colorbar
    c = colorbar;
    c.Label.String = 'Epoch Number';
    clim([min(epoch_subset), max(epoch_subset)]);
    colormap(ax, turbo);
    
    % Set axis properties
    grid on;
    box on;
    axis equal;
    
    % Format axes with appropriate scaling
    ax.XAxis.Exponent = 7;
    ax.YAxis.Exponent = 7;
    
    xlabel('X (m)');
    ylabel('Y (m)');
    title('X-Y View of Satellite Constellation');
    
    hold off;
end

function plotChiSquareStatistics(chi_square, criticalValue, numEpochs, ax)
    % Plot chi-square test statistics over time
    
    % Create axes
    axes(ax);
    hold on;
    
    % Generate x-axis (epoch numbers)
    epochs = 1:numEpochs;
    
    % Plot chi-square values
    plot(epochs, chi_square, 'b-', 'LineWidth', 1.5);
    
    % Add critical value line
    yline(criticalValue, 'r--', 'Critical Value', 'LineWidth', 1.5);
    
    % Highlight regions above critical value
    idx_above = chi_square > criticalValue;
    if any(idx_above)
        area(epochs(idx_above), chi_square(idx_above), 'FaceColor', 'r', ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    
    % Set axis properties
    grid on;
    box on;
    xlabel('Epoch');
    ylabel('Chi-Square Value');
    title('RAIM Fault Detection Test Statistics');
    
    hold off;
end

function plotProtectionLevels(protectionLevels, numEpochs, ax)
    % Plot protection levels over time
    
    % Create axes
    axes(ax);
    hold on;
    
    % Generate x-axis (epoch numbers)
    epochs = 1:numEpochs;
    
    % Plot protection levels
    plot(epochs, protectionLevels, 'g-', 'LineWidth', 2);
    
    % Fill area under the curve
    area(epochs, protectionLevels, 'FaceColor', [0.4, 0.8, 0.4], ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Set axis properties
    grid on;
    box on;
    xlabel('Epoch');
    ylabel('Protection Level (m)');
    title('RAIM Protection Levels');
    ylim([0, max(protectionLevels)*1.1]);
    
    hold off;
end