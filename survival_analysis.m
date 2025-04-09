clear;
close all;
load Time_passed_Seizure_all
load Dur_seizure_all

% Initialize empty arrays to hold the combined data
T = [];           % Survival times (time to next seizure)
E = [];           % Event indicators (1 = seizure occurred)
Duration = [];    % Seizure durations
TimeBefore = [];  % Time before seizures
RatID = [];       % Identifier for each rat

for rat = 1:6
    % Initialize temporary variables to hold data for each rat
    combined_time_before = [];
    combined_seizure_durations = [];
    
    for day = 1:40
        % Extract time before seizures and durations for the current rat and day
        time_before_seizures = Time_passed_Seizure_all{rat, day};
        seizure_durations = Dur_seizure_all{rat, day};
        
        % Concatenate the data for the current day
        combined_time_before = [combined_time_before; time_before_seizures];
        combined_seizure_durations = [combined_seizure_durations; seizure_durations];
    end
    
    num_seizures = length(combined_time_before);
    
    for s = 1:num_seizures - 1  % Exclude the last seizure if the next seizure time is unknown
        % Survival time is the time to next seizure
        survival_time =  combined_time_before(s+1);
        T = [T; survival_time];
        
        % Event indicator (1 = seizure occurred)
        E = [E; 1];
        
        % Covariate: duration of current seizure
        Duration = [Duration; combined_seizure_durations(s)];
        
        % Covariate: time before the current seizure
        TimeBefore = [TimeBefore; combined_time_before(s)];
        
        % Identifier for the rat
        RatID = [RatID; rat];
    end
end

% Create a table for the Cox model
data = table(T, E, Duration, TimeBefore, RatID);

% Fit the Cox model separately for each rat and store the hazard ratios
HR_rat = zeros(6, 2);  % Store hazard ratios for each rat (Duration and TimeBefore)
CI_lower = zeros(6, 2);
CI_upper = zeros(6, 2);

for rat = 1:6
    % Filter data for each rat
    rat_idx = (RatID == rat);
    
    % Fit the Cox model with both covariates (Duration and TimeBefore)
    [b, ~, ~, stats] = coxphfit([data.Duration(rat_idx), data.TimeBefore(rat_idx)], data.T(rat_idx), ...
                                'Censoring', ~data.E(rat_idx)); % Note: This line still uses the E vector
    
    % Extract hazard ratios
    HR_rat(rat, :) = exp(b);  % Exponentiated coefficients give the hazard ratios
    
    % Confidence intervals (CI)
    CI_lower(rat, :) = exp(b - 1.96 * stats.se);
    CI_upper(rat, :) = exp(b + 1.96 * stats.se);
end

% Plot hazard ratios for each rat (with confidence intervals)
figure;
for i = 1:2  % i = 1 for Duration, i = 2 for TimeBefore
    subplot(1, 2, i);
    errorbar(1:6, HR_rat(:, i), HR_rat(:, i) - CI_lower(:, i), CI_upper(:, i) - HR_rat(:, i), 'o', 'LineWidth', 2);
    ylim([0.9 1.1])
    xticks(1:6);
    xticklabels({'Rat T2', 'Rat T3', 'Rat T4', 'Rat T5', 'Rat T6', 'Rat T7'});
    
    if i == 1
        title('Hazard Ratios for Seizure Duration');
    else
        title('Hazard Ratios for Time Before Seizure');
    end
    
    xlabel('Rat');
    ylabel('Hazard Ratio');
    grid on;
end

% Kaplan-Meier Survival Curve without censoring
figure;
hold on;
for rat = 1:6
    % Filter data for each rat
    rat_idx = (RatID == rat);
    
    % Compute the ECDF (survivor function) for each rat
    % No need for censoring argument since we assume all events are observed
    [ecdf_y, ecdf_x] = ecdf(T(rat_idx), 'Function', 'survivor');
    
    % Plot the survival curve for each rat
    stairs(ecdf_x, ecdf_y, 'LineWidth', 2, 'DisplayName', sprintf('Rat T%d', rat+1));
end
xlabel('Time to Next Seizure');
ylabel('Survival Probability');
title('Kaplan-Meier Survival Curve by Rat');
legend show;
grid on;
