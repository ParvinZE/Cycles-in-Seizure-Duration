% Load data for on_periods and off_periods
clear
close all
all_rat_seizure_counts_on_periods=load('all_rat_seizure_counts_on_periods.mat'); % Load matrix for seizure counts on_periods
all_rat_seizure_counts_off_periods=load('all_rat_seizure_counts_off_periods.mat'); % Load matrix for seizure counts off_periods
 all_rat_seizure_durations_on_periods=load('all_rat_seizure_durations_on_periods.mat'); % Load cell for seizure durations on_periods
all_rat_seizure_durations_off_periods=load('all_rat_seizure_durations_off_periods.mat'); % Load cell for seizure durations off_periods

% Perform comparison for seizure counts
% Compute summary statistics
mean_counts_on = mean(all_rat_seizure_counts_on_periods.all_rat_seizure_counts, 1);
mean_counts_off = mean(all_rat_seizure_counts_off_periods.all_rat_seizure_counts, 1);
std_counts_on = std(all_rat_seizure_counts_on_periods.all_rat_seizure_counts, 0, 1);
std_counts_off = std(all_rat_seizure_counts_off_periods.all_rat_seizure_counts, 0, 1);

% Perform statistical test
[h_counts, p_counts] = ttest2(all_rat_seizure_counts_on_periods.all_rat_seizure_counts, all_rat_seizure_counts_off_periods.all_rat_seizure_counts, 'dim', 1);



% Perform comparison for seizure durations
% Compute summary statistics
mean_durations_on = cellfun(@mean, all_rat_seizure_durations_on_periods.all_rat_seizure_durations);
mean_durations_off = cellfun(@mean, all_rat_seizure_durations_off_periods.all_rat_seizure_durations);
std_durations_on = cellfun(@std, all_rat_seizure_durations_on_periods.all_rat_seizure_durations);
std_durations_off = cellfun(@std, all_rat_seizure_durations_off_periods.all_rat_seizure_durations);

% Perform statistical test
[h_durations, p_durations] = ttest2(cellfun(@mean, all_rat_seizure_durations_off_periods.all_rat_seizure_durations),cellfun(@mean, all_rat_seizure_durations_on_periods.all_rat_seizure_durations), 'dim', 1);


% Define custom colors for the boxplots
color_on = [0.2 0.6 0.2]; % Green for on_periods
color_off = [0.8 0.2 0.2]; % Red for off_periods



mean_durations_off_mean=nanmean(mean_durations_off,1);
mean_durations_on_mean=nanmean(mean_durations_on,1);



% Define the number of rats
num_rats = 6;



% Adjust figure properties
sgtitle('Comparison of Seizure Durations between On Periods and Off Periods');


% Define quantiles (e.g., lowest and highest)
quantile_low = 0.25; % Lowest quantile
quantile_high = 0.75; % Highest quantile

% Create a figure for the histograms
figure;

% Iterate over each rat
for rat_idx = 1:num_rats
    % Extract seizure durations for on_periods and off_periods for the current rat
    durations_on = cell2mat(all_rat_seizure_durations_on_periods.all_rat_seizure_durations(:, rat_idx));
    durations_off = cell2mat(all_rat_seizure_durations_off_periods.all_rat_seizure_durations(:, rat_idx));
    
    % Compute quantiles for on_periods and off_periods
    quantile_on_low = quantile(durations_on, quantile_low);
    quantile_on_high = quantile(durations_on, quantile_high);
    quantile_off_low = quantile(durations_off, quantile_low);
    quantile_off_high = quantile(durations_off, quantile_high);
    
    % Create subplots for the current rat
    subplot(num_rats, 1, rat_idx);
    
    % Plot histograms for on_periods and off_periods
    histogram(durations_on, 'BinWidth', 5, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;
    histogram(durations_off, 'BinWidth', 5, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    % Plot dashed lines for quantiles
    plot([quantile_on_low, quantile_on_low], ylim, '--b');
    plot([quantile_on_high, quantile_on_high], ylim, '--b');
    plot([quantile_off_low, quantile_off_low], ylim, '--r');
    plot([quantile_off_high, quantile_off_high], ylim, '--r');


%     
    
    % Set plot title
    title(['Rat T' num2str(rat_idx+1)]);
    
    % Set axis labels
    xlabel('Duration of Seizures (seconds)');
    ylabel('Frequency');
    
    % Add legend
    legend('On Periods', 'Off Periods', 'Location', 'best');
end

% Adjust figure properties
sgtitle('Comparison of Seizure Durations between On Periods and Off Periods');


% Define custom colors for the boxes
color_on = [0.2 0.2 0.8];   % Blue color
color_off = [0.8 0.2 0.2];  % Red color

% Define custom colors for the boxes
color_on = [0.2 0.2 0.8];   % Blue color
color_off = [0.8 0.2 0.2];  % Red color


% Load seizure duration data for on periods
load('all_rat_seizure_durations_on_periods.mat'); % Ensure this loads the correct structure
durations_on = all_rat_seizure_durations_on_periods.all_rat_seizure_durations;

% Load seizure duration data for off periods
load('all_rat_seizure_durations_off_periods.mat'); % Ensure this loads the correct structure
durations_off = all_rat_seizure_durations_off_periods.all_rat_seizure_durations;

% Define colors similar to seaborn's husl palette
colors = [
    0.193, 0.678, 0.211; % Greenish
    0.867, 0.520, 0.227; % Orangeish
    0.121, 0.308, 0.757; % Blueish
    0.749, 0.225, 0.407; % Reddish
    0.184, 0.184, 0.671; % Darker blue
    0.988, 0.755, 0.164  % Yellowish
];

% Number of rats
num_rats = size(durations_on, 2);

% Create a figure for the histograms
% Load seizure duration data for on periods
load('all_rat_seizure_durations_on_periods.mat'); % Ensure this loads the correct structure
durations_on = all_rat_seizure_durations_on_periods.all_rat_seizure_durations;

% Load seizure duration data for off periods
load('all_rat_seizure_durations_off_periods.mat'); % Ensure this loads the correct structure
durations_off = all_rat_seizure_durations_off_periods.all_rat_seizure_durations;

% Define a blue color
dark_blue_color = [0, 0, 0.5]; % RGB values for dark blue

% Number of rats
num_rats = size(durations_on, 2);

% Create a figure for the histograms
figure;

% Iterate over each rat
for rat_idx = 1:num_rats
    % Extract seizure durations for on_periods and off_periods for the current rat
    durations_on_rat = durations_on(:, rat_idx);
    durations_off_rat = durations_off(:, rat_idx);
    
    % Combine durations for on and off periods
    combined_durations = [cell2mat(durations_on_rat); cell2mat(durations_off_rat)];
    
    % Create a subplot for the current rat
    subplot(num_rats, 1, rat_idx);
    
    % Plot histogram for combined durations
    histogram(combined_durations, 'BinWidth', 3, 'FaceColor', dark_blue_color, ...
              'EdgeColor', 'white', 'FaceAlpha',0.8);
    
    % Set plot title and axis labels
    title(['Rat T' num2str(rat_idx + 1)]);
    xlabel('Duration of Seizures (seconds)');
    ylabel('Frequency');
    
    % Add legend
   % legend('Durations');
    
    % Perform statistical test comparing on and off durations (optional)
    [h_durations, p_durations,stats_durations] = ttest2(cell2mat(durations_on_rat), cell2mat(durations_off_rat), 'Vartype', 'unequal');
     t_value = stats_durations(1);
     df=stats_durations(2);
    disp(['Rat T' num2str(rat_idx + 1) ' p-value: ' num2str(p_durations)]);
    fprintf('t(%0.1f) = %0.2f, p = %0.4f\n', df, t_value, p_durations);
end

% Adjust figure properties
%sgtitle('Comparison of Combined Seizure Durations for On and Off Periods');

