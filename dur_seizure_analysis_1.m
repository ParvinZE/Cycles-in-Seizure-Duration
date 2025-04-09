clear
close all

% Initialize variables for 6 rats across 40 days
Dur_seizure_all = cell(6, 40);
Num_seizures_all = zeros(6, 40);
Time_passed_Seizure_all = cell(6, 40);

% Initialize mean/median duration containers
mean_dur_all = zeros(6, 40);
median_dur_all = zeros(6, 40);

Dur_seizure = [];
count_short = 0;
count_long = 0;
Short_seizures = {};
Long_seizures = {};
base_path = 'data'; % or './data', or 'path/to/your/dataset'

for rat = 2:7
    for day = 1:40
        duration_file = fullfile(base_path, sprintf('T%d', rat), ...
            sprintf('duration_T%d_ictal_Day%d.mat', rat, day));
        timepassed_file = fullfile(base_path, sprintf('T%d', rat), ...
            sprintf('TimePassed_T%d_ictal_Day%d.mat', rat, day));

        load(duration_file);
        load(timepassed_file);

        if isempty(Duration)
            Duration = nan;
        end

        Dur_seizure_all{rat-1, day} = Duration;
        Time_passed_Seizure_all{rat-1, day} = cell2mat(TimePassed);
        Num_seizures_all(rat-1, day) = length(Duration);

        if ~isempty(Duration)
            mean_dur_all(rat-1, day) = nanmean(Duration);
            median_dur_all(rat-1, day) = nanmedian(Duration);
        end
    end
end

% Calculate mean and median seizure counts across rats
mean_num_seizures = nanmean(Num_seizures_all, 1);
median_num_seizures = nanmedian(Num_seizures_all, 1);

% Calculate mean and median seizure durations across rats
mean_dur_seizures = zeros(1, 40);
median_dur_seizures = zeros(1, 40);
for day = 1:40
    all_durations = [];
    for rat = 2:7
        all_durations = [all_durations; Dur_seizure_all{rat-1, day}];
    end
    mean_dur_seizures(day) = mean(all_durations, 'omitnan');
    median_dur_seizures(day) = median(all_durations, 'omitnan');
end

% Compute SEM for number and duration of seizures
std_error_num_seizures = nanstd(Num_seizures_all, 1) / sqrt(size(Num_seizures_all, 1));
std_error_dur_seizures = zeros(1, 40);
for day = 1:40
    all_durations = [];
    for rat = 2:7
        all_durations = [all_durations; Dur_seizure_all{rat-1, day}];
    end
    std_error_dur_seizures(day) = nanstd(all_durations) / sqrt(sum(~isnan(all_durations)));
end

% ==== PLOT: Mean & SEM for Number and Duration of Seizures ====
figure;

% Plot SEM for number of seizures (shaded)
std_error_fill_num = [mean_num_seizures + std_error_num_seizures; ...
                      mean_num_seizures - std_error_num_seizures];
fill([1:40, fliplr(1:40)], [std_error_fill_num(1,:), fliplr(std_error_fill_num(2,:))], ...
     [0.9 0.9 0.9], 'EdgeColor', 'none');
hold on;

% Mean number of seizures
yyaxis left;
plot(1:40, mean_num_seizures, 'k-', 'LineWidth', 2, 'DisplayName', 'Mean Number of Seizures');
ylabel('Number of Seizures');

% SEM for duration of seizures (shaded)
std_error_fill_dur_upper = [];
std_error_fill_dur_lower = [];
for day = 1:40
    std_error = std_error_dur_seizures(day);
    mean_dur = mean_dur_seizures(day);
    if ~isnan(std_error) && ~isnan(mean_dur)
        std_error_fill_dur_upper(day) = mean_dur + std_error;
        std_error_fill_dur_lower(day) = mean_dur - std_error;
    end
end

yyaxis right;
ax = gca;
ax.YAxis(2).Color = 'b';
fill([1:40, fliplr(1:40)], ...
     [std_error_fill_dur_upper, fliplr(std_error_fill_dur_lower)], ...
     [0.3 0.3 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on;

% Mean duration of seizures
plot(1:40, mean_dur_seizures, 'b-', 'LineWidth', 2, 'DisplayName', 'Mean Duration of Seizures');
ylabel('Duration of Seizures (seconds)');
xlabel('Day');
title('Mean of Number and Duration of Seizures Per Day');
xlim([4, 40]);
legend({'SEM of Number of Seizures','Mean Number of Seizures', ...
        'SEM of Duration of Seizures', 'Mean Duration of Seizures'}, ...
        'Location', 'best');
grid on;


% % Create box plot



% === More detailed grouping based on time before and after ===
durations_long_before = [];
durations_long_until = [];
durations_long_before_short = [];
durations_short_until = [];

for rat = 2:7
    durations_all = [];
    timepassed_all = [];

    for day = 1:40
        if isempty(Time_passed_Seizure_all{rat-1, day})
            Time_passed_Seizure_all{rat-1, day} = nan;
        end

        durations = Dur_seizure_all{rat-1, day};
        timepassed = Time_passed_Seizure_all{rat-1, day};

        if ~any(isnan(timepassed)) && ~any(isnan(durations))
            durations_all = [durations_all; durations(:)];
            timepassed_all = [timepassed_all; timepassed(:)];
        end
    end

    timepassed_all_after = [timepassed_all(2:end); NaN];

    durations_long_before = [durations_long_before; durations_all(timepassed_all >= 1 * 24 * 60 * 60)];
    durations_long_before_short = [durations_long_before_short; durations_all(timepassed_all < 1.5 * 60)];
    durations_long_until = [durations_long_until; durations_all(timepassed_all_after(1:end-1) >= 1 * 24 * 60 * 60)];
    durations_short_until = [durations_short_until; durations_all(timepassed_all_after(1:end-1) < 1.5 * 60)];
end

% === Boxplot customization ===
color_before = [0.1, 0.6, 0.1];        % Green
color_until = [0.1, 0.1, 0.8];         % Blue
color_before_short = [1, 0.5, 0];      % Orange
color_until_short = [1, 0, 0];         % Red

combined_data = [durations_long_before; durations_long_until; ...
                 durations_long_before_short; durations_short_until];

combined_labels = [repmat("Time Before Seizure >= 1 day", length(durations_long_before), 1); ...
                   repmat("Time Until Next Seizure >= 1 day", length(durations_long_until), 1); ...
                   repmat("Time Before Seizure < 1.5 minutes", length(durations_long_before_short), 1); ...
                   repmat("Time Until Next Seizure < 1.5 minutes", length(durations_short_until), 1)];

colors = {color_before, color_until, color_before_short, color_until_short};

figure;
boxplot_handle = boxplot(combined_data, combined_labels, 'Whisker', Inf);
ylabel('Duration (seconds)');
grid on;

% Fill box colors
h = findobj(gca, 'Tag', 'Box');
for k = 1:length(h)
    patch(get(h(k), 'XData'), get(h(k), 'YData'), colors{k}, ...
          'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Style lines
median_handles = findobj(gca, 'Tag', 'Median');
whisker_handles = findobj(gca, 'Tag', 'Whisker');
cap_handles = findobj(gca, 'Tag', 'Cap');

for k = 1:length(colors)
    if k <= length(median_handles)
        set(median_handles(k), 'Color', colors{k}, 'LineWidth', 1.5);
    end
    if 2*k-1 <= length(whisker_handles)
        set(whisker_handles(2*k-1:2*k), 'Color', colors{k}, 'LineWidth', 1.5);
    end
    if 2*k-1 <= length(cap_handles)
        set(cap_handles(2*k-1:2*k), 'Color', colors{k}, 'LineWidth', 1.5);
    end
end

% Outlier color
outliers = findobj(gca, 'Tag', 'Outliers');
set(outliers, 'MarkerEdgeColor', 'k');

% Add n labels
num_points = [length(durations_long_before), length(durations_long_until), ...
              length(durations_long_before_short), length(durations_short_until)];
x_positions = 1:length(num_points);
for k = 1:length(num_points)
    text(x_positions(k), max(combined_data) - 0.05 * range(combined_data), ...
         ['n = ' num2str(num_points(k))], ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');
end

% Storage
delay_before_i_long_then_next_delay_long = [];
delay_before_i_long_then_next_delay_short = [];

for rat = 2:7
    durations_all = [];
    timepassed_all = [];

    for day = 1:40
        durations = Dur_seizure_all{rat-1, day};
        timepassed = Time_passed_Seizure_all{rat-1, day};
        if isempty(durations) || isempty(timepassed)
            continue
        end
        durations_all = [durations_all; durations(:)];
        timepassed_all = [timepassed_all; timepassed(:)];
    end

    for i = 1:length(timepassed_all)  % go only up to i+1
        % If seizure i happened after a long delay
        if timepassed_all(i) >= 1*24*60*60
            % Look at the next seizure (i+1), and check its delay BEFORE it
            delay = timepassed_all(i);  % delay before seizure i+1
            if delay >= 1*24*60*60
                delay_before_i_long_then_next_delay_long = [delay_before_i_long_then_next_delay_long; timepassed_all(i+1)];
                end
        end
    end
end


% Combine the two sets of durations
combined_durations = [ delay_before_i_long_then_next_delay_long];

% Create group labels
group_labels = [repmat("Next Seizure Preceded by ≥1 day", length(delay_before_i_long_then_next_delay_long), 1)];

% Plot boxplot

% Define custom colors
color_next = [0.5, 0.1, 0.5]; % Purple-like
color_before  = [0.1, 0.1, 0.8]; % Blue-like
box_colors = { color_next,color_before};

% Combine data
combined_delays = [delay_before_i_long_then_next_delay_long/3600];
group_labels = [repmat("Inter-seizure interval until next seizure preceded by ≥1 day", length(delay_before_i_long_then_next_delay_long), 1)];

% Create boxplot
figure;
boxplot(combined_delays, group_labels, 'Whisker', Inf);
ylabel('Inter-Seizure intervals (hours)');
%title('Seizures Preceded by Long Delay and the Next Seizure');
grid on;

% Color customization using patch
h = findobj(gca, 'Tag', 'Box');
for k = 1:length(h)
    xData = get(h(k), 'XData');
    yData = get(h(k), 'YData');
    patch(xData, yData, box_colors{k}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Style medians, whiskers, caps
median_handles = findobj(gca, 'Tag', 'Median');
whisker_handles = findobj(gca, 'Tag', 'Whisker');
cap_handles = findobj(gca, 'Tag', 'Cap');

for k = 1:length(box_colors)
    if k <= length(median_handles)
        set(median_handles(k), 'Color', box_colors{k}, 'LineWidth', 1.5);
    end
    if 2*k-1 <= length(whisker_handles)
        set(whisker_handles(2*k-1:2*k), 'Color', box_colors{k}, 'LineWidth', 1.5);
    end
    if 2*k-1 <= length(cap_handles)
        set(cap_handles(2*k-1:2*k), 'Color', box_colors{k}, 'LineWidth', 1.5);
    end
end

% Style outliers
outliers = findobj(gca, 'Tag', 'Outliers');
set(outliers, 'MarkerEdgeColor', 'k');

% Add n = count above each box
num_points = [length(delay_before_i_long_then_next_delay_long)];
x_positions = 1:length(num_points);
for k = 1:length(num_points)
    text(x_positions(k), max(combined_delays) - 0.05 * range(combined_durations), ...
         ['n = ' num2str(num_points(k))], ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');
end

