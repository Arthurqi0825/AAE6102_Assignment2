function mask = combined_skymask_gnss()
    % combined_skymask_gnss - Combines skymask visualization with GNSS positioning.
    %--------------------------------------------------------------------------
    % 1. USER INPUT / LOAD DATA
    %--------------------------------------------------------------------------
    % Ground-truth approx position (lat, lon, alt)
    lat0 = 22.3198722;
    lon0 = 114.209101777778;
    alt0 = 3.0;
    
    % Define relative paths (files are in the same directory as this script)
    navSolutionFile = 'navSolutionResults.mat';
    skymaskFile = 'skymask_A1_urban.csv';
    
    % Load GNSS data
    if exist(navSolutionFile, 'file')
        navSolutions1 = load(navSolutionFile);
        mask.pr = navSolutions1.navSolutions.correctedP;
        mask.az = navSolutions1.navSolutions.az;
        mask.el = navSolutions1.navSolutions.el;
        mask.pos = navSolutions1.navSolutions.satPositions;
    else
        error('File not found: %s', navSolutionFile);
    end
    disp('Loaded Data Sizes:');
    disp(size(mask.pr));
    disp(size(mask.az));
    disp(size(mask.el));
    disp(size(mask.pos));
    
    % Load skymask CSV
    M = readmatrix(skymaskFile);
    mask.maskAz = M(:,1);
    mask.maskEl = M(:,2);
    maskElVec = nan(361,1);
    for i = 1:numel(mask.maskAz)
        a = round(mod(mask.maskAz(i),360));
        maskElVec(a+1) = mask.maskEl(i);
    end
    idx = find(~isnan(maskElVec));
    mask.maskElVec = interp1(idx-1, maskElVec(idx), (0:360)','linear','extrap');
    
    % Plot skymask
    figure;
    subplot(1, 2, 1); % 1 row, 2 columns, first subplot
    azimuth_rad = deg2rad(mask.maskAz); % Convert azimuth degrees to radians
    polarplot(azimuth_rad, 90 - mask.maskEl, 'LineWidth', 2, 'Color', 'b');
    title('Skymask Visualization (Polar)');
    rlim([0 90]); % Limit radius to 90 degrees
    thetalim([0 360]); % Cover full azimuth range
    grid on;
    
    subplot(1, 2, 2); % 1 row, 2 columns, second subplot
    plot(mask.maskAz, mask.maskEl, 'LineWidth', 2, 'Color', 'b');
    title('Skymask Visualization (Horizontal)');
    xlabel('Azimuth Angle (degrees)');
    ylabel('Elevation Angle (degrees)');
    ylim([0 90]); % Limit elevation to 90 degrees
    grid on;
    sgtitle('Skymask Visualizations'); % Overall title for both subplots
    
    % WGS-84 ellipsoid
    mask.wgs = wgs84Ellipsoid('meters');
    [x0, y0, z0] = geodetic2ecef(mask.wgs, lat0, lon0, alt0);
    dt0 = 0;
    c = 299792458;
    [nSat, nEpoch] = size(mask.pr);
    sol = nan(nEpoch, 4);
    nVis = zeros(nEpoch, 1);
    
    %--------------------------------------------------------------------------
    % 2. MAIN LOOP OVER EPOCHS
    %--------------------------------------------------------------------------
    for k = 1:nEpoch
        mask.rho_k = mask.pr(:, k);
        mask.az_k = mask.az(:, k);
        mask.el_k = mask.el(:, k);
        mask.Psat = squeeze(mask.pos(:, :, k))';
        
        % Skymask visibility & weights
        vis = false(nSat, 1);
        w = zeros(nSat, 1);
        for i = 1:nSat
            a = mod(mask.az_k(i), 360);
            ai = floor(a) + 1;
            el_block = mask.maskElVec(ai);
            if mask.el_k(i) > el_block
                w(i) = sin(deg2rad(mask.el_k(i) - el_block)) * sin(deg2rad(mask.el_k(i)));
                vis(i) = true;
            end
        end
        idx = find(vis);
        nVis(k) = numel(idx);
        fprintf('Visible satellites at epoch %d: %d\n', k, nVis(k));
        if numel(idx) < 4
            continue; % cannot solve
        end
        
        % WLS via Gauss–Newton
        x_est = [x0; y0; z0; dt0];
        tol = 1e-4;
        for iter = 1:10
            m = numel(idx);
            mask.H = zeros(m, 4);
            mask.r = zeros(m, 1);
            mask.W = diag(w(idx));
            for ii = 1:m
                i = idx(ii);
                rho_hat = norm(mask.Psat(i, :)' - x_est(1:3));
                pred = rho_hat + c * x_est(4);
                mask.r(ii) = mask.rho_k(i) - pred;
                u = (x_est(1:3) - mask.Psat(i, :)') / rho_hat;
                mask.H(ii, 1:3) = u';
                mask.H(ii, 4) = -c;
            end
            % Solve for Δx
            dx = (mask.H' * mask.W * mask.H) \ (mask.H' * mask.W * mask.r);
            x_est = x_est + dx;
            if norm(dx) < tol, break; end
        end
        sol(k, :) = x_est';
        x0 = x_est(1);
        y0 = x_est(2);
        z0 = x_est(3);
        dt0 = x_est(4);
    end
    
    %--------------------------------------------------------------------------
    % 3. POST-PROCESS & PLOT
    %--------------------------------------------------------------------------
    fprintf('Visible sats per epoch (min/max): %d / %d\n', min(nVis), max(nVis));
    valid = find(~any(isnan(sol), 2));
    if isempty(valid)
        disp('No valid solutions to convert to geodetic coordinates.');
        return;
    end
    [lat, lon, ~] = ecef2geodetic(mask.wgs, sol(valid, 1), sol(valid, 2), sol(valid, 3));
    if isempty(lat) || isempty(lon)
        disp('No valid latitude and longitude values to plot.');
        return;
    end
    
    % Plot results
    figure;
    plot(lon, lat, 'r.', 'MarkerSize', 10); % Estimated positions
    hold on;
    plot(lon0, lat0, 'bx', 'MarkerSize', 10); % Ground truth
    xlabel('Longitude (°)');
    ylabel('Latitude (°)');
    legend('Estimated Position', 'Ground Truth');
    title('GNSS Positioning in Urban Environment');
    grid on;
end