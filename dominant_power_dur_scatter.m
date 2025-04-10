clear
close all

% Define the full path to the rat's folder
rat_data_folder =  'data';

% Initialize variables to store power and durations for all rats
all_rat_seizure_durations = [];
all_rat_dominant_powers = [];
all_rat_dominant_frequencies = [];
rat_names = {};
rat_colors = lines(6); % Predefined color scheme for 6 rats

% Initialize an array to store correlation coefficients for each rat
rat_correlations = [];
p_values = [];

% Define logarithmic frequency bands up to 130 Hz
fs = 2048; % Sampling frequency
max_freq = 130; % Maximum frequency of interest (not the Nyquist frequency, but the limit for our analysis)
num_bands = 10; % Number of frequency bands
frequency_bands = logspace(log10(1), log10(max_freq), num_bands); % Logarithmically spaced bands up to 130 Hz
band_names = strcat(string(round(frequency_bands(1:end-1))), '-', string(round(frequency_bands(2:end))), ' Hz'); % Create labels for bands

% Filter rat_folders to include only directories representing rat folders
rat_folders = dir(rat_data_folder);
rat_folders = rat_folders([rat_folders.isdir]); % Keep only directories
is_rat_folder = arrayfun(@(x) startsWith(x.name, 'T') && ~strcmp(x.name, '.') && ~strcmp(x.name, '..'), rat_folders);
rat_folders = rat_folders(is_rat_folder);

% Figure for individual plots
figure('Name', 'Dominant Frequency Band vs. Duration for Each Rat');
tiledlayout(3, 2); % 3x2 grid for subplots

% Prepare figure for combined plot
figure('Name', 'Dominant Frequency Band vs. Duration (All Rats)');
hold on;
rat_Dominant_Powers = {};

% Iterate through each rat's folder
for file2_idx = 1:numel(rat_folders)
    rat_folder = rat_folders(file2_idx).name;
    rat_name = rat_folder; % Get the rat name
    rat_names{end+1} = rat_name; % Store the rat name for coloring in the combined plot
    rat_idx = file2_idx;
    
    disp(['Processing ' rat_name]);
    
    % Define the full path to the rat's folder
    rat_folder_path = fullfile(rat_data_folder, rat_folder);
    
    % Initialize variables to store data for this rat
    rat_seizure_durations = [];
    rat_dominant_powers = [];
    rat_dominant_frequencies = [];
    
    % Iterate over each day in the rat's folder
    day_files = dir(fullfile(rat_folder_path, '*.mat')); % Assuming all files are .mat files
    
    for file_idx = 1:numel(day_files)
        day_file = day_files(file_idx).name;
        if startsWith(day_file, ['EEG_' rat_name '_ictal_Day']) % Matching the EEG files
            % Load EEG signal
            eeg_data = load(fullfile(rat_folder_path, day_file));
            eeg_signal = eeg_data.Seizure_signal; % Assuming the EEG signal is stored under this variable
            
            % Convert cell array to numeric array for each seizure
            for seizure_idx = 1:length(eeg_signal)
                seizure = eeg_signal{seizure_idx};
                
                if isempty(seizure)
                    continue; % Skip empty EEG signals
                end
                
                % Calculate the duration of the seizure using the length of the EEG signal
                seizure_duration = length(seizure) / fs; % Adjusted to the sampling rate
                
                % Initialize variables to store power in each logarithmic band
                band_powers = zeros(length(frequency_bands)-1, 1);
                
                % Calculate power for each logarithmic frequency band (up to 130 Hz)
                for band_idx = 1:length(frequency_bands) - 1
                    f_low = frequency_bands(band_idx);
                    f_high = min(frequency_bands(band_idx + 1), max_freq); % Ensure f_high does not exceed 130 Hz
                    
                    % Skip invalid frequency ranges (e.g., if f_low >= f_high)
                    if f_low >= f_high
                        continue;
                    end
                    
                    % Calculate power in the current frequency range using Welch's method (manual calculation of PSD)
                    [pxx, f] = pwelch(seizure, [], [], [], fs); % Power spectral density using Welch's method
                    band_power_idx = find(f >= f_low & f <= f_high); % Find the frequency indices within the band
                    band_powers(band_idx) = sum(pxx(band_power_idx)); % Sum the power in this frequency range
                end
                
                % Find the dominant frequency band (max power)
                [dominant_power, dominant_band_idx] = max(band_powers);
                dominant_frequency_band = band_names{dominant_band_idx}; % Get the name of the dominant band
                
                % Store the seizure duration and dominant power/frequency for this rat
                rat_seizure_durations = [rat_seizure_durations; seizure_duration];
                rat_dominant_powers = [rat_dominant_powers; dominant_power];
                rat_dominant_frequencies = [rat_dominant_frequencies; dominant_band_idx]; % Store index for color coding
                
            end
        end
    end
    
    % Append rat-specific data to the overall data
    all_rat_seizure_durations = [all_rat_seizure_durations; rat_seizure_durations];
    all_rat_dominant_powers = [all_rat_dominant_powers; rat_dominant_powers];
    all_rat_dominant_frequencies = [all_rat_dominant_frequencies; rat_dominant_frequencies];
    rat_Dominant_Powers{rat_idx} = rat_dominant_powers;
    
    % Calculate the correlation for this rat
    if ~isempty(rat_seizure_durations) && ~isempty(rat_dominant_powers)
        [r, p] = corr(rat_seizure_durations, rat_dominant_powers, 'Type', 'Pearson');
        rat_correlations = [rat_correlations; r];
        p_values = [p_values; p];
    end
    
    % Plot for this specific rat in the subplot
    figure(1);
    nexttile;
    scatter(rat_seizure_durations, rat_dominant_powers, 'filled', 'MarkerFaceColor', rat_colors(rat_idx, :));
    xlabel('Seizure Duration (s)');
    ylabel('Dominant EEG Power');
    title(['Rat ' rat_name]);
    
    % Plot data for this rat in the combined plot with different color
    figure(2);
    scatter(rat_seizure_durations, rat_dominant_powers, 'filled', 'MarkerFaceColor', rat_colors(rat_idx, :));
end

% Final touches on combined plot
figure(2);
xlabel('Seizure Duration (s)');
ylabel('Dominant EEG Power');
title('Relationship between Seizure Duration and Dominant EEG Power (All Rats)');
legend(rat_names, 'Location', 'best');
hold off;

% Apply Bonferroni correction to p-values if applicable
if ~isempty(p_values)
    corrected_p_values = p_values * numel(p_values); % Bonferroni correction
    
    % Ensure the correlation, p-values, and corrected p-values are column vectors
    if size(rat_correlations, 2) > 1
        rat_correlations = rat_correlations';
    end
    if size(p_values, 2) > 1
        p_values = p_values';
    end
    if size(corrected_p_values, 2) > 1
        corrected_p_values = corrected_p_values';
    end
    
    % Check if the dimensions match
    if length(rat_correlations) == length(p_values) && length(p_values) == length(corrected_p_values)
        % Proceed to create the table only if sizes match
        results_table = array2table([rat_correlations, p_values, corrected_p_values], ...
            'VariableNames', {'Correlation_Coefficient', 'P_value', 'Corrected_P_value'}, ...
            'RowNames', rat_names);
        
        % Convert table to a figure
        figure('Name', 'Correlation and P-values');
        uitable('Data', table2cell(results_table), ...
                'ColumnName', results_table.Properties.VariableNames, ...
                'RowName', results_table.Properties.RowNames, ...
                'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
        
        % Save the figure
        exportgraphics(gcf, 'Correlation_Pvalues_Table.png', 'Resolution', 300);
    else
        error('The dimensions of rat_correlations, p_values, and corrected_p_values do not match.');
    end
else
    disp('No p-values were calculated. Skipping the p-value table creation.');
end


