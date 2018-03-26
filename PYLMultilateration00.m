function [GEO_position]=PYLMultilateration00
    %% Multilateration Looped Random Error Analysis Simulation
    % Creator Name: Jared Burbey
    % Date Created: 2016-07-15
    % Last Date Modified: 2016-09-16
    % Modified By: Griffin Curtis
    % Modifications: Adapted to work with TDOA simulation, fixed
    % 3d plot errors including est coordinates were centered around 0 error (should
    % have been the mean of the error), and the 3d plot only showed one
    % standard deviation in each direction, not 3.
    %
    % References: Dr. Roggemann, Ph.D in Electro-Optics, EERC 503 MTU
    % campus
    %
    % This function calculates the position of a Geostationary satellite
    % using multilateration with 3 ground receivers and 1 LEO space
    % receiver.
    % All distance mearsurements are in meters unless stated otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Initialize for fmincon
    clc;
    close all
    global AvgLoc
    global Geo_Positions
    global est_tdoa_with_noise1
    global est_tdoa_with_noise2
    global est_tdoa_with_noise3
    global est_tdoa_with_noise4
    global est_tdoa_with_noise5
    global est_tdoa_with_noise6
    global R_a_lla
    global R_b_lla
    global R_c_lla
    global R_d_lla
    global R_a 
    global R_b 
    global R_c
    global R_d
    global GEO_coordinates
    global GEO_lla
    
    iterations = 1000;
    calculated_position = zeros(iterations, 3);
    errors = zeros(iterations,3);
    direct_errors = zeros(1,iterations);

    % The horizontal RMS value from the OEM615 datasheet
    horiz_err_gnd = 1.5;
    horiz_err_leo = 10; % 10 Modify this to test different conditions for error budget
    
    vert_horiz_err_ratio = 1.9; % Ranges from about 1.4 (optimistic) to 1.9 (pessimistic)
    
    vert_err_gnd = horiz_err_gnd * vert_horiz_err_ratio;
    vert_err_leo = horiz_err_leo * vert_horiz_err_ratio;
    
    horiz_err_gnd = horiz_err_gnd / sqrt(2);
    horiz_err_leo = horiz_err_leo / sqrt(2);
    
    % The std dev error in each direction for ground and LEO
    err_gnd = [horiz_err_gnd horiz_err_gnd vert_err_gnd];
    err_leo = [horiz_err_leo horiz_err_leo vert_err_leo];
    
    for i = 1:1:iterations
    
        % Geostationary satellite position and inital position guess
        GEO_true_position = GEO_coordinates;
        GEO_guess_lla = mean([R_a_lla; R_b_lla; R_c_lla; R_d_lla]);
        GEO_guess_lla(3) = GEO_lla(3) + 1000000*randn();
        initial_guess = lla2ecef(GEO_guess_lla);
        
        % timing error of receivers
        errbudget = 22;  %22 90
        t_error = (errbudget*randn(1,4))/1e9;
        
        err_a = (err_leo .* randn(1,3)) * lla2ecefneu(R_a_lla);
        err_b = (err_gnd .* randn(1,3)) * lla2ecefneu(R_b_lla);
        err_c = (err_gnd .* randn(1,3)) * lla2ecefneu(R_c_lla);
        err_d = (err_gnd .* randn(1,3)) * lla2ecefneu(R_d_lla);
        
        % position error of receiver. Each row is a receiver. Each column
        % is a dimension in ecef
        p_error = [err_a; err_b; err_c; err_d];

        % fmincon options that set max interations, function tolerances, and
        % removes a debug display output
        options = optimoptions('lsqnonlin', 'TolX', 1e-20, 'TolFun', 1e-20, 'MaxIter', 1e4, 'MaxFunEvals',1e4,'Jacobian', 'on', 'Display','none', 'Algorithm', 'levenberg-marquardt', 'ScaleProblem', 'jacobian');

%% Begin fmincon and Error Calculations
        % Calaculate the optimized solution that returns the position for the
        % GEO Sat
        GEO_position = lsqnonlin(@MLAT,initial_guess,[],[],options);

        % Calculate the error of the fmincon output
        error = GEO_position - GEO_true_position;
        direct_error = sqrt(error(1)^2 + error(2)^2 + error(3)^2);
        direct_error_km = direct_error/1000; % error in meters turned into kilometers
        
        errors(i,1:3) = error/1000; % each column is the error in x,y,z respectively
        direct_errors(i) = direct_error_km;
%% Radical New Code for Beam Mapping
        Geo_Positions(i,1:3) = GEO_position; 
        AvgLocX = mean(Geo_Positions(:,1));
        AvgLocY = mean(Geo_Positions(:,2));
        AvgLocZ = mean(Geo_Positions(:,3));
        AvgLoc = [AvgLocX, AvgLocY, AvgLocZ];
        AvgLoc = ecef2lla(AvgLoc,'WGS84');
        
        % Create Standard Deviation Matrix for x,y,z
        known_position = GEO_true_position/1000; % in km
        calculated_position(i,1:3) = GEO_position/1000; % in km

       end
    
    %% General data statistics
    % Plot the data after multiple iterations
    figure(1);
    histogram(direct_errors, 'Normalization', 'probability');
    title('Frequency of occurence for distance error');
    xlabel('Direct Distance Error, km');
    ylabel('Frequency');
    
    % Output average data
    sorted_errors = sort(direct_errors);
    fprintf('\nMean = %.4f km\nMedian = %.4f km\n\n',mean(direct_errors),median(sorted_errors));
    
    %% Plot the success sphere and error volume
    % Create 3D Ellipsoid of the Failure condition sphere
    figure(2);
    [x_true,y_true,z_true] = ellipsoid(0,0,0,50,50,50);
    surf(x_true,y_true,z_true,'FaceAlpha',.2,'EdgeAlpha',.1);
    axis equal
    hold on
    
    %%
    % Create 3D ellipsoid for the actual errors (3-sigma)
    covar = cov(errors); % This is the covariance matrix. We need to make this transform from std-deviation space to error-space
    center_pos = mean(errors);
    [x_err,y_err,z_err] = ellipsoid(0,0,0,3,3,3, 8); % The x,y,z positions of the points on the surface of the ellipsoid. Right now it's a 3-sigma ellipsoid
    
    % To take the covariance (sigma^2) and turn it into standard deviation,
    % we have to do an eigen-decomposition and square-root the eigenvalues
    [covar_V, covar_D] = eig(covar);
    covar_D = sqrt(covar_D); % Take the sqare root of the elements in the diagonal eigenvalue matrix
    covar = covar_V * covar_D / covar_V; % covar is now a std-deviation matrix instead of a covariance matrix!
    
    % Now we transform the x,y,z from the ellipsoid from std-deviation
    % space to error space
    for i = 1:numel(x_err)
        v = covar * [x_err(i); y_err(i); z_err(i)] + transpose(center_pos);
        x_err(i) = v(1);
        y_err(i) = v(2);
        z_err(i) = v(3);
    end
    
    surf(x_err,y_err,z_err,'LineStyle','none','FaceColor','r');
    
    %% Magic over, carry on
    % Create scatter plot of actual data
    figure(3)
    hold on
    surf(x_true,y_true,z_true,'FaceAlpha',.2,'EdgeAlpha',.1);
    axis equal
    scatter3(errors(1:iterations,1),errors(1:iterations,2),errors(1:iterations,3), '.');
    hold off
    title('Calculated Volume vs Maximum Allowed Volume of the GEO satellite''s position');
    xlabel('Distance X, km');
    ylabel('Distance Y, km');
    zlabel('Distance Z, km');
    
    
    figure(4)
    hold on
    surf(x_err,y_err,z_err,'FaceAlpha',0.3);
    scatter3(errors(1:iterations,1),errors(1:iterations,2),errors(1:iterations,3), '.');
    xlabel('Distance X, km');
    ylabel('Distance Y, km');
    zlabel('Distance Z, km');
    hold off
    
    % Plot the bell curve of the results
    mu = mean(errors) * covar_V;
    sigma = diag(covar_D);
    
    [~, maj_index] = max(sigma); % index of the maximum standard deviation. Major axis
    [~, min_index] = min(sigma); % index of the minimum standard deviation. Minor axis
    int_index = 6 - (min_index + maj_index); % index of the middle standard deviation. Intermediate axis
    
    mu = mu([maj_index int_index min_index]);
    sigma = sigma([maj_index int_index min_index]);
    
    range_bellcurve = -35:.1:35;
    norm_1 = normpdf(range_bellcurve,mu(1),sigma(1));
    norm_2 = normpdf(range_bellcurve,mu(2),sigma(2));
    norm_3 = normpdf(range_bellcurve,mu(3),sigma(3));

    figure
    plot(range_bellcurve,norm_1);
    title('Major axis error distribution');
    xlabel('Distance, km');
    ylabel('Probability of occurence');
    
    figure
    plot(range_bellcurve,norm_2);
    title('Intermediate axis error distribution');
    xlabel('Distance, km');
    ylabel('Probability of occurence');
    
    figure
    plot(range_bellcurve,norm_3);
    title('Minor axis error distribution');
    xlabel('Distance, km');
    ylabel('Probability of occurence');
    
    % Save positions to a .mat file
    save('multilateration_positions.mat','known_position','calculated_position');
        
    fprintf('Major axis        mean, std. dev = %+.4f, %.4f', mu(1), sigma(1));
    fprintf('\nIntermediate axis mean, std. dev = %+.4f, %.4f', mu(2), sigma(2));
    fprintf('\nMinor axis        mean, std. dev = %+.4f, %.4f\n', mu(3), sigma(3));
    NotUsedVar=multbeammap(AvgLoc, GEO_lla);
%% Multilateration error of a point p
    % This function calculates the error and jacobian matrix for the system
    % of equations that specify the position of the GEOsat.
    % The minimum point of this function is the location of the GEOsat
    function [F,J] = MLAT(p)
        c = 299792458; % m/s speed of light
        
        % Time differences in m
        t_ba = (est_tdoa_with_noise1 + t_error(1) + t_error(2)) * c;
        t_ca = (est_tdoa_with_noise2 + t_error(1) + t_error(3)) * c;
        t_da = (est_tdoa_with_noise3 + t_error(1) + t_error(4)) * c;
        t_cb = (est_tdoa_with_noise4 + t_error(2) + t_error(3)) * c;
        t_db = (est_tdoa_with_noise5 + t_error(2) + t_error(4)) * c;
        t_dc = (est_tdoa_with_noise6 + t_error(3) + t_error(4)) * c;

        % Positions in m
        p_a = R_a + p_error(1,1:3);
        p_b = R_b + p_error(2,1:3);
        p_c = R_c + p_error(3,1:3);
        p_d = R_d + p_error(4,1:3);
        
        % Difference between the receiver positions and the GEOsat position
        d_a = p_a - p;
        d_b = p_b - p;
        d_c = p_c - p;
        d_d = p_d - p;
        
        % Distances between receivers and the GEOsat Position
        n_a = norm(d_a);
        n_b = norm(d_b);
        n_c = norm(d_c);
        n_d = norm(d_d);
        
        % The system of equations. Ideally, all should equal 0 and be
        % minimized at the position of the GEOsat
        F = [n_b - n_a - t_ba;
             n_c - n_a - t_ca;
             n_d - n_a - t_da;
             n_c - n_b - t_cb;
             n_d - n_b - t_db;
             n_d - n_c - t_dc];
         
        if nargout > 1
            % The matrix of partial derivatives of the above. Speeds up the
            % minimization routine
            J = [d_a(1)/n_a - d_b(1)/n_b, d_a(2)/n_a - d_b(2)/n_b, d_a(3)/n_a - d_b(3)/n_b;
                 d_a(1)/n_a - d_c(1)/n_c, d_a(2)/n_a - d_c(2)/n_c, d_a(3)/n_a - d_c(3)/n_c;
                 d_a(1)/n_a - d_d(1)/n_d, d_a(2)/n_a - d_d(2)/n_d, d_a(3)/n_a - d_d(3)/n_d;
                 d_b(1)/n_b - d_c(1)/n_c, d_b(2)/n_b - d_c(2)/n_c, d_b(3)/n_b - d_c(3)/n_c;
                 d_b(1)/n_b - d_d(1)/n_d, d_b(2)/n_b - d_d(2)/n_d, d_b(3)/n_b - d_d(3)/n_d;
                 d_c(1)/n_c - d_d(1)/n_d, d_c(2)/n_c - d_d(2)/n_d, d_c(3)/n_c - d_d(3)/n_d];
        end
    end

%% Calculates the north, east, up matrix
    % This function takes a position in lla and calculates the directions
    % of north, east, and up and returns them as normalized row vectors in 
    % a matrix

    function neu = lla2ecefneu(lla)
        dl = 0.001;
        da = 1;
        
        % Finite differences
        n = lla2ecef(lla + [dl, 0, 0]) - lla2ecef(lla);
        e = lla2ecef(lla + [0, dl, 0]) - lla2ecef(lla);
        u = lla2ecef(lla + [0, 0, da]) - lla2ecef(lla);
        
        % Normalize
        n = n / norm(n);
        e = e / norm(e);
        u = u / norm(u);
        
        % Make the matrix
        neu = [n; e; u];
    end
end