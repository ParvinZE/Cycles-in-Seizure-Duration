clear
close all

% Define the full path to the rat's folder
rat_data_folder = 'data';

% Initialize variables to store power, durations, and dominant frequencies for all rats
all_rat_dominant_frequencies = [];
rat_labels = [];
rat_names = {};
rat_colors = lines(6); % Predefined color scheme for 6 rats

% Define frequency bands up to 130 Hz
fs = 2048; % Sampling frequency
max_freq = 130; % Maximum frequency of interest (not the Nyquist frequency)
num_bands =10; % Number of frequency bands
frequency_bands = logspace(log10(1), log10(max_freq), num_bands); % Logarithmic bands up to 130 Hz

% Filter to include only rat folders
rat_folders = dir(rat_data_folder);
rat_folders = rat_folders([rat_folders.isdir]); % Keep only directories
is_rat_folder = arrayfun(@(x) startsWith(x.name, 'T') && ~strcmp(x.name, '.') && ~strcmp(x.name, '..'), rat_folders);
rat_folders = rat_folders(is_rat_folder);

% Process each rat folder
rat_dominant_frequencies = cell(numel(rat_folders), 1); % To store dominant frequencies for each rat
rat_power_dominant_frequencies=cell(numel(rat_folders), 1);
for rat_idx = 1:numel(rat_folders)
    rat_folder = rat_folders(rat_idx).name;
    rat_names{end+1} = rat_folder; % Store the rat name
    disp(['Processing ' rat_folder]);

    % Path to the rat's folder
    rat_folder_path = fullfile(rat_data_folder, rat_folder);

    % Initialize variables to store data for this rat
    rat_dominant_frequencies{rat_idx} = [];

    % Iterate over each day in the rat's folder
    day_files = dir(fullfile(rat_folder_path, '*.mat')); % Assuming all files are .mat files
    
    for file_idx = 1:numel(day_files)
        day_file = day_files(file_idx).name;
        if startsWith(day_file, ['EEG_' rat_folder '_ictal_Day']) % Matching the EEG files
            % Load EEG signal
            eeg_data = load(fullfile(rat_folder_path, day_file));
            eeg_signal = eeg_data.Seizure_signal; % Assuming the EEG signal is stored here
            
            % Convert cell array to numeric array for each seizure
            for seizure_idx = 1:length(eeg_signal)
                seizure = eeg_signal{seizure_idx};
                
                if isempty(seizure)
                    continue; % Skip empty EEG signals
                end
                
                % Initialize variables to store power in each logarithmic band
                band_powers = zeros(length(frequency_bands)-1, 1);

                % Calculate power for each logarithmic frequency band (up to 130 Hz)
                [pxx, f] = pwelch(seizure, [], [], [], fs); % Power spectral density
                for band_idx = 1:length(frequency_bands) - 1
                    f_low = frequency_bands(band_idx);
                    f_high = min(frequency_bands(band_idx + 1), max_freq);
                    
                    if f_low >= f_high
                        continue; % Skip invalid frequency ranges
                    end
                    
                    band_power_idx = find(f >= f_low & f <= f_high);
                    band_powers(band_idx) = sum(pxx(band_power_idx));
                end
                
                % Find the dominant frequency band (max power)
                [max_power, dominant_band_idx] = max(band_powers);
                dominant_frequency_band = frequency_bands(dominant_band_idx); % Get dominant frequency

                % Save dominant frequency for this seizure
                rat_dominant_frequencies{rat_idx} = [rat_dominant_frequencies{rat_idx}; dominant_frequency_band];
                rat_power_dominant_frequencies{rat_idx} = [rat_power_dominant_frequencies{rat_idx}; max_power];
            end
        end
    end
    
    % Append the data and create group labels for boxplot
    all_rat_dominant_frequencies = [all_rat_dominant_frequencies; rat_dominant_frequencies{rat_idx}];
    rat_labels = [rat_labels; repmat(rat_names(rat_idx), numel(rat_dominant_frequencies{rat_idx}), 1)];
end

% Create a single figure for the combined box plot
figure('Name', 'Dominant Frequency Box Plots for All Rats', 'Position', [100, 100, 1200, 600]);

% Create the grouped boxplot for all rats
boxplot(all_rat_dominant_frequencies, rat_labels, 'Colors', lines(numel(rat_names)), 'Symbol', 'o', ...
        'BoxStyle', 'outline', 'Widths', 0.7, 'Notch', 'on', 'Whisker', 1.5);

% Set labels and title
xlabel('Rats');
ylabel('Dominant Frequency (Hz)');
title('Dominant Frequency Distribution for All Rats');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

% Color the boxplot using patch
boxes = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxes)
    patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), rat_colors(j, :), 'FaceAlpha', 0.3);
end

% Add grid for better readability
grid on;
