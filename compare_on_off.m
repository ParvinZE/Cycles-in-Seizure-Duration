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

% % Create boxplots
% figure;
% subplot(2, 1, 1);
% boxplot(all_rat_seizure_counts_on_periods.all_rat_seizure_counts, 'Labels', {'Rat T2', 'Rat T3', 'Rat T4', 'Rat T5', 'Rat T6', 'Rat T7'});
% title('Seizure Counts - On Periods');
% 
% subplot(2, 1, 2);
% boxplot(all_rat_seizure_counts_off_periods.all_rat_seizure_counts, 'Labels',{'Rat T2', 'Rat T3', 'Rat T4', 'Rat T5', 'Rat T6', 'Rat T7'});
% title('Seizure Counts - Off Periods');

% Perform comparison for seizure durations
% Compute summary statistics
mean_durations_on = cellfun(@mean, all_rat_seizure_durations_on_periods.all_rat_seizure_durations);
mean_durations_off = cellfun(@mean, all_rat_seizure_durations_off_periods.all_rat_seizure_durations);
std_durations_on = cellfun(@std, all_rat_seizure_durations_on_periods.all_rat_seizure_durations);
std_durations_off = cellfun(@std, all_rat_seizure_durations_off_periods.all_rat_seizure_durations);

% Perform statistical test
[h_durations, p_durations] = ttest2(cellfun(@mean, all_rat_seizure_durations_off_periods.all_rat_seizure_durations),cellfun(@mean, all_rat_seizure_durations_on_periods.all_rat_seizure_durations), 'dim', 1);

% % Create boxplots
% figure;
% subplot(2, 1, 1);
% boxplot(cellfun(@mean, all_rat_seizure_durations_off_periods.all_rat_seizure_durations), 'Labels', {'Rat 1', 'Rat 2', 'Rat 3', 'Rat 4', 'Rat 5', 'Rat 6'});
% title('Seizure Durations - On Periods');
% 
% subplot(2, 1, 2);
% boxplot(cellfun(@mean, all_rat_seizure_durations_on_periods.all_rat_seizure_durations), 'Labels', {'Rat 1', 'Rat 2', 'Rat 3', 'Rat 4', 'Rat 5', 'Rat 6'});
% title('Seizure Durations - Off Periods');
% 
% % Display statistical results
% disp('Statistical Results for Seizure Counts:');
% disp(['p-value for seizure counts: ' num2str(p_counts)]);
% disp(['Hypothesis test for seizure counts (H = 0: No significant difference, H = 1: Significant difference): ' num2str(h_counts)]);
% disp(' ');
% disp('Statistical Results for Seizure Durations:');
% disp(['p-value for seizure durations: ' num2str(p_durations)]);
% disp(['Hypothesis test for seizure durations (H = 0: No significant difference, H = 1: Significant difference): ' num2str(h_durations)]);
% 

% Define custom colors for the boxplots
color_on = [0.2 0.6 0.2]; % Green for on_periods
color_off = [0.8 0.2 0.2]; % Red for off_periods

% % Create figure for seizure counts
% figure;
% 
% % Plot boxplots for seizure counts
% 
% for rat_idx = 1:6
%     x_on = rat_idx - 0.2; % Shift for on_periods
%     x_off = rat_idx + 0.2; % Shift for off_periods
% 
%     % Plot boxplots for on_periods and off_periods side by side
%     boxplot([all_rat_seizure_counts_on_periods.all_rat_seizure_counts(:, rat_idx), all_rat_seizure_counts_off_periods.all_rat_seizure_counts(:, rat_idx)], 'Colors', [color_on; color_off], 'Widths', 0.4, 'Positions', [x_on, x_off]);
%     hold on;
% end
% title('Seizure Counts');
% xlabel('Rat');
% ylabel('Counts');
% xticks(1:6);
% xticklabels({'Rat T2', 'Rat T3', 'Rat T4', 'Rat T5', 'Rat T6', 'Rat T7'});
% legend('On Periods', 'Off Periods', 'Location', 'best');
% hold off;

% % Create figure for seizure durations
% figure;
% 
% % Plot boxplots for seizure durations
% 
% for rat_idx = 1:6
%     x_on = rat_idx - 0.2; % Shift for on_periods
%     x_off = rat_idx + 0.2; % Shift for off_periods
% 
%     % Plot boxplots for on_periods and off_periods side by side
%     boxplot([cellfun(@mean, all_rat_seizure_durations_on_periods.all_rat_seizure_durations(:, rat_idx)), cellfun(@mean, all_rat_seizure_durations_off_periods.all_rat_seizure_durations(:, rat_idx))], 'Colors', [color_on; color_off],'Widths', 0.4, 'Positions', [x_on, x_off]);
%     hold on;
% end
% title('Seizure Durations');
% xlabel('Rat');
% ylabel('Durations');
% xticks(1:6);
% xticklabels({'Rat T2', 'Rat T3', 'Rat T4', 'Rat T5', 'Rat T6', 'Rat T7'});
% legend('On Periods', 'Off Periods', 'Location', 'best');
% hold off;
% 
% 

mean_durations_off_mean=nanmean(mean_durations_off,1);
mean_durations_on_mean=nanmean(mean_durations_on,1);



% Define the number of rats
num_rats = 6;

% % Create a figure for the histograms
% figure;
% 
% % Iterate over each rat
% for rat_idx = 1:num_rats
%     % Extract seizure durations for on_periods and off_periods for the current rat
%     durations_on = all_rat_seizure_durations_on_periods.all_rat_seizure_durations(:, rat_idx);
%     durations_off = all_rat_seizure_durations_off_periods.all_rat_seizure_durations(:, rat_idx);
% 
%     % Create subplots for the current rat
%     subplot(num_rats, 1, rat_idx);
% 
%     % Plot histograms for on_periods and off_periods
%     histogram(cell2mat(durations_on), 'BinWidth', 5, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     hold on;
%     histogram(cell2mat(durations_off), 'BinWidth', 5, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% 
%     % Set plot title
%     title(['Rat T' num2str(rat_idx+1)]);
% 
%     % Set axis labels
%     xlabel('Duration of Seizures (seconds)');
%     ylabel('Frequency');
% 
%     % Add legend
%     legend('On Periods', 'Off Periods');
%     % Perform statistical test
% [h_durations, p_durations] = ttest2(cell2mat(durations_on),cell2mat(durations_off), 'dim', 1);
% disp(p_durations)
% end

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

%      % Add text for quantile values
%     text(quantile_on_low, 0.9*max(ylim), ['Q' num2str(quantile_low) ' = ' num2str(quantile_on_low)], 'Color', 'b');
%     text(quantile_on_high, 0.9*max(ylim), ['Q' num2str(quantile_high) ' = ' num2str(quantile_on_high)], 'Color', 'b');
%     text(quantile_off_low, 0.8*max(ylim), ['Q' num2str(quantile_low) ' = ' num2str(quantile_off_low)], 'Color', 'r');
%     text(quantile_off_high, 0.8*max(ylim), ['Q' num2str(quantile_high) ' = ' num2str(quantile_off_high)], 'Color', 'r');
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

% figure;
% 
% % Define custom whisker lengths (adjust as needed)
% whisker_length = 1.5;
% 
% % Define the spacing between boxes for different rats
% rat_spacing = 1.5; % Adjust as needed
% 
% % Plot boxplots for seizure durations
% for rat_idx = 1:6
%     % Adjust the x positions for each rat
%     x_on = rat_idx * rat_spacing - 0.2; % Shift for on_periods
%     x_off = rat_idx * rat_spacing + 0.2; % Shift for off_periods
% 
%     % Extract seizure durations for on_periods and off_periods for the current rat
%     durations_on = cell2mat(all_rat_seizure_durations_on_periods.all_rat_seizure_durations(:, rat_idx));
%     durations_off = cell2mat(all_rat_seizure_durations_off_periods.all_rat_seizure_durations(:, rat_idx));
% 
%     % Plot boxplots for on_periods and off_periods side by side
%     bp1 = boxplot(durations_on, 'Colors', 'b', 'Widths', 0.6, 'Notch', 'on', 'Positions', x_on);
%     hold on;
%     bp2 = boxplot(durations_off, 'Colors', 'r', 'Widths', 0.6, 'Notch', 'on', 'Positions', x_off);
%     hold on;
% 
%     % Fill the inside of the boxes with colors
%     h1 = findobj(bp1, 'Tag', 'Box');
%     h2 = findobj(bp2, 'Tag', 'Box');
%     patch(get(h1, 'XData'), get(h1, 'YData'), 'b', 'FaceAlpha', 0.5); % Fill on_periods box
%     patch(get(h2, 'XData'), get(h2, 'YData'), 'r', 'FaceAlpha', 0.5); % Fill off_periods box
% 
%     % Get handles to the outliers
%     outliers_on = findobj(bp1, 'Tag', 'Outliers');
%     outliers_off = findobj(bp2, 'Tag', 'Outliers');
% 
%     % Change the color of the outliers
%     set(outliers_on, 'MarkerEdgeColor', 'b');
%     set(outliers_off, 'MarkerEdgeColor', 'r');
% end
% 
% % Customize the plot
% title('Seizure Durations');
% xlabel('Rat');
% ylabel('Durations');
% legend('On Periods', 'Off Periods', 'Location', 'best');
% xticks((1+0.5):rat_spacing:(6+0.5)*rat_spacing); % Set x-axis ticks to the middle of each group of boxes
% xticklabels({'Rat T2', 'Rat T3', 'Rat T4', 'Rat T5', 'Rat T6', 'Rat T7'});
% grid on;
% hold off;


% figure;
% 
% % Define the spacing between boxes for different rats
% rat_spacing = 1.5; % Adjust as needed
% 
% % Custom box height scaling factor (use Inf to extend boxes to cover all data points)
% box_scale_factor = Inf;
% 
% % Plot boxplots for seizure durations
% for rat_idx = 1:6
%     % Adjust the x positions for each rat
%     x_on = rat_idx * rat_spacing - 0.2; % Shift for on_periods
%     x_off = rat_idx * rat_spacing + 0.2; % Shift for off_periods
% 
%     % Extract seizure durations for on_periods and off_periods for the current rat
%     durations_on = cell2mat(all_rat_seizure_durations_on_periods.all_rat_seizure_durations(:, rat_idx));
%     durations_off = cell2mat(all_rat_seizure_durations_off_periods.all_rat_seizure_durations(:, rat_idx));
% 
%     % Check if durations are not empty for on_periods and off_periods
%     if ~isempty(durations_on)
%         % Calculate quartiles and IQR for on periods
%         q1_on = quantile(durations_on, 0.25);
%         q3_on = quantile(durations_on, 0.75);
%         iqr_on = q3_on - q1_on;
% 
%         % Adjust box height to include all data points
%         box_low_on = min(durations_on);  % Set lower limit of box to the minimum value
%         box_high_on = max(durations_on); % Set upper limit of box to the maximum value
% 
%         % Manually plot box rectangle (blue for on_periods)
%         fill([x_on-0.3 x_on+0.3 x_on+0.3 x_on-0.3], [box_low_on box_low_on box_high_on box_high_on], 'b', 'FaceAlpha', 0.5);
% 
%         % Plot the medians as horizontal lines inside the boxes
%         median_on = median(durations_on);
%         plot([x_on-0.3 x_on+0.3], [median_on median_on], 'k-', 'LineWidth', 2);
%         hold on
%         % Plot individual data points for on periods
%         scatter(repmat(x_on, size(durations_on)), durations_on, 'b', 'filled');
%         hold on
%     end
% 
%     if ~isempty(durations_off)
%         % Calculate quartiles and IQR for off periods
%         q1_off = quantile(durations_off, 0.25);
%         q3_off = quantile(durations_off, 0.75);
%         iqr_off = q3_off - q1_off;
% 
%         % Adjust box height to include all data points
%         box_low_off = min(durations_off);  % Set lower limit of box to the minimum value
%         box_high_off = max(durations_off); % Set upper limit of box to the maximum value
% 
%         % Manually plot box rectangle (red for off_periods)
%         fill([x_off-0.3 x_off+0.3 x_off+0.3 x_off-0.3], [box_low_off box_low_off box_high_off box_high_off], 'r', 'FaceAlpha', 0.5);
% 
%         % Plot the medians as horizontal lines inside the boxes
%         median_off = median(durations_off);
%         plot([x_off-0.3 x_off+0.3], [median_off median_off], 'k-', 'LineWidth', 2);
%         hold on
%         % Plot individual data points for off periods
%         scatter(repmat(x_off, size(durations_off)), durations_off, 'r', 'filled');
%     end
%     hold on
% end
% 
% % Customize the plot
% title('Seizure Durations');
% xlabel('Rat');
% ylabel('Durations');
% legend('On Periods', 'Off Periods', 'Location', 'best');
% xticks((1+0.5):rat_spacing:(6+0.5)*rat_spacing); % Set x-axis ticks to the middle of each group of boxes
% xticklabels({'Rat T2', 'Rat T3', 'Rat T4', 'Rat T5', 'Rat T6', 'Rat T7'});
% grid on;
% hold off;
% 


% figure;
% % Extract seizure durations for on_periods and off_periods for the current rat
% for rat_idx=1:6
%     % Compute mean and SEM for on_periods and off_periods
%     mean_on = cellfun(@mean, all_rat_seizure_durations_on_periods.all_rat_seizure_durations(:, rat_idx));
%     mean_off =cellfun(@mean, all_rat_seizure_durations_off_periods.all_rat_seizure_durations(:, rat_idx));
%     std_on = cellfun(@std, all_rat_seizure_durations_on_periods.all_rat_seizure_durations(:, rat_idx));
%     std_off =cellfun(@std, all_rat_seizure_durations_off_periods.all_rat_seizure_durations(:, rat_idx));
% 
%     sem_on = std_on / sqrt(size(durations_on, 1)); % Standard error of the mean for on_periods
% 
%     sem_off = std_off / sqrt(size(durations_off, 1)); % Standard error of the mean for off_periods
% 
% 
% 
% 
%     % Create subplot for the current rat
%     subplot(6, 1, rat_idx);
% 
%     % Plot mean and SEM for on_periods
%     errorbar(1:length(mean_on), mean_on, sem_on, 'ob', 'LineWidth', 1.5);
%     hold on;
% 
%     % Plot mean and SEM for off_periods
%     errorbar(1:length(mean_off), mean_off, sem_off, 'or', 'LineWidth', 1.5);
% 
%     % Customize subplot
%     xlabel('Day');
%     ylabel('Seizure Duration');
%     title(['Rat T' num2str(rat_idx+1)]);
%     legend('On Periods', 'Off Periods');
%     grid on;
% end

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
    [h_durations, p_durations] = ttest2(cell2mat(durations_on_rat), cell2mat(durations_off_rat), 'Vartype', 'unequal');
    disp(['Rat T' num2str(rat_idx + 1) ' p-value: ' num2str(p_durations)]);
end

% Adjust figure properties
%sgtitle('Comparison of Combined Seizure Durations for On and Off Periods');
