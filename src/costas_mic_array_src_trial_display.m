function costas_mic_array_src_trial_display(errors, frequencies, snr_values, placements, position_sets, source_positions, all_predictions, output_dir)
    % This function calculates and displays statistics from the errors data,
    % and visualizes microphone placements and predictions.

    % Ensure output directory exists
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    columns_placement = length(placements);

    %% Average errors for each placement, SNR, and frequency
    manual_errors_snr_freq = zeros(columns_placement, length(snr_values) * length(frequencies));
    column_names_snr_freq = cell(size(manual_errors_snr_freq, 2), 1);
    col_idx = 1;
    for snr_idx = 1:length(snr_values)
        for freq_idx = 1:length(frequencies)
            for placement_idx = 1:columns_placement
                manual_errors_snr_freq(placement_idx, col_idx) = mean(errors(placement_idx, :, snr_idx, freq_idx), 'all');
            end
            column_names_snr_freq{col_idx} = sprintf('SNR_%ddB_Freq_%dHz', snr_values(snr_idx), frequencies(freq_idx));
            col_idx = col_idx + 1;
        end
    end

    T_snr_freq = array2table(manual_errors_snr_freq, 'RowNames', placements, 'VariableNames', column_names_snr_freq);
    T_snr_freq = highlight_min_with_star(T_snr_freq);

    %% Average errors across all frequencies for each SNR per placement
    manual_errors_snr = zeros(columns_placement, length(snr_values));
    column_names_snr = cell(length(snr_values), 1);
    for snr_idx = 1:length(snr_values)
        for placement_idx = 1:columns_placement
            manual_errors_snr(placement_idx, snr_idx) = mean(errors(placement_idx, :, snr_idx, :), [2, 4]);
        end
        column_names_snr{snr_idx} = sprintf('SNR_%ddB', snr_values(snr_idx));
    end
    T_snr = array2table(manual_errors_snr, 'RowNames', placements, 'VariableNames', column_names_snr);
    T_snr = highlight_min_with_star(T_snr);

    %% Average errors across all SNRs for each frequency per placement
    manual_errors_freq = zeros(columns_placement, length(frequencies));
    column_names_freq = cell(length(frequencies), 1);
    for freq_idx = 1:length(frequencies)
        for placement_idx = 1:columns_placement
            manual_errors_freq(placement_idx, freq_idx) = mean(errors(placement_idx, :, :, freq_idx), [2, 3]);
        end
        column_names_freq{freq_idx} = sprintf('Freq_%dHz', frequencies(freq_idx));
    end
    T_freq = array2table(manual_errors_freq, 'RowNames', placements, 'VariableNames', column_names_freq);
    T_freq = highlight_min_with_star(T_freq);

    %% Overall stats (mean, median, variance, and 75th percentile) for each placement
    manual_stats = zeros(columns_placement, 5);
    stat_names = {'Mean', 'Median', 'Variance', '75th Percentile','RMSE'};
    for placement_idx = 1:columns_placement
        placement_errors = errors(placement_idx, :, :, :);
        placement_errors = placement_errors(:); % Flatten
        manual_stats(placement_idx, 1) = mean(placement_errors);
        manual_stats(placement_idx, 2) = median(placement_errors);
        manual_stats(placement_idx, 3) = var(placement_errors);
        manual_stats(placement_idx, 4) = prctile(placement_errors, 75);
        %Focus on Large Errors: Squaring emphasizes larger errors, making RMSE sensitive to significant deviations
        manual_stats(placement_idx, 5) = sqrt(mean(placement_errors.^2)); % Compute RMSE

    end
    T_stats = array2table(manual_stats, 'RowNames', placements, 'VariableNames', stat_names);
    T_stats = highlight_min_with_star(T_stats);

    %% Display results
    disp('Error Table: Per Placement, Per SNR x Frequency Pair:');
    disp(T_snr_freq);

    disp('Error Table: Per Placement, Averaged Across Frequencies:');
    disp(T_freq);

    disp('Error Table: Per Placement, Averaged Across SNRs:');
    disp(T_snr);

    disp('Overall Error Stats (Mean, Median, Variance, 75th Percentile):');
    disp(T_stats);

    %% Box plot of errors
    all_errors = [];
    group_labels = [];
    for placement_idx = 1:columns_placement
        placement_errors = errors(placement_idx, :, :, :);
        placement_errors = placement_errors(:); % Flatten
        all_errors = [all_errors; placement_errors];
        group_labels = [group_labels; repmat(placements(placement_idx), length(placement_errors), 1)];
    end

    figure;
    boxplot(all_errors, group_labels, 'Labels', placements);
    hold on;

    % Add average markers for each placement
    unique_labels = unique(group_labels);
    placement_avgs = zeros(1, length(unique_labels));
    for i = 1:length(unique_labels)
        placement_avg = mean(all_errors(strcmp(group_labels, unique_labels{i})));
        placement_avgs(i) = placement_avg;
        plot(i, placement_avg, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
        text(i, placement_avg, sprintf('%.2f', placement_avg), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'blue');
    end

    % Highlight the smallest average error
    [~, min_idx] = min(placement_avgs);
    plot([min_idx - 0.4, min_idx + 0.4], [placement_avgs(min_idx), placement_avgs(min_idx)], 'g-', 'LineWidth', 2);

    % Add legend placeholders (single entry for Average and Minimum Average)
    h_avg = plot(nan, nan, 'r*', 'MarkerSize', 10, 'DisplayName', 'Average');
    h_min_avg = plot(nan, nan, 'g-', 'LineWidth', 2, 'DisplayName', 'Minimum Average');

    legend([h_avg, h_min_avg], {'Average', 'Minimum Average'}, 'Location', 'best');

    %% 3D Visualization of Microphone Arrangements, Sources, and Predictions
    figure;
hold on;
%markers = {'o', '^', 's', 'd', 'p', 'h', 'x', '+', 'v', '<'};
%colors = {'b', 'm', 'g', 'c', 'y', 'k', 'r', 'b', 'm', 'g'};
    markers = {'o', '^', 's', 'd', 'p', 'h', 'x', '+', 'v','<',};
%    colors = {'b', 'm', 'g', 'c', 'y', 'k', 'r', 'b', 'm', 'g'};
colors = {'b', 'm', 'g', 'c', 'y', 'k', 'r', ...  % Standard colors
          [0.5, 0.2, 0.8], [0.9, 0.3, 0.3], [0.3, 0.7, 0.2], ... % Custom RGB colors
          [0.1, 0.5, 0.9], [0.6, 0.2, 0.5], [0.9, 0.7, 0.1], ...
          [0.3, 0.3, 0.3], [0.7, 0.4, 0.1], [0.2, 0.8, 0.8]};
legend_handles = [];
legend_labels = {};

for p = 1:length(position_sets)
    array_positions = position_sets{p};
    h = scatter3(array_positions(:, 1), array_positions(:, 2), array_positions(:, 3), 100, ...
                 colors{mod(p-1, length(colors))+1}, markers{mod(p-1, length(markers))+1}, ... %'filled', ...
                     'MarkerFaceColor', colors{mod(p-1, length(colors))+1}, ... % Ensure face color matches
                     'MarkerEdgeColor', 'k'); % Explicitly define edge color
                 
    legend_handles = [legend_handles, h];
    legend_labels{end+1} = placements{p};

    for s = 1:size(source_positions, 1)
        source_position = source_positions(s, :);
        for t_idx = 1:size(all_predictions, 2)
            predicted_position = all_predictions{p, t_idx, s};
            scatter3(predicted_position(1), predicted_position(2), predicted_position(3), 50, ...
                     colors{mod(p-1, length(colors))+1}, markers{mod(p-1, length(markers))+1}, ...
                     'MarkerFaceColor', colors{mod(p-1, length(colors))+1}, ... % Ensure face color matches
                     'MarkerEdgeColor', 'k'); % Explicitly define edge color
        end
    end
end

%legend(legend_handles, legend_labels);
%hold off;

%{
    figure;
    hold on;
    markers = {'o', '^', 's', 'd', 'p', 'h', 'x', '+', 'v','<',};
    colors = {'b', 'm', 'g', 'c', 'y', 'k', 'r', 'b', 'm', 'g'};

    legend_handles = [];
    legend_labels = {};

    for p = 1:length(position_sets)
        array_positions = position_sets{p};
        h = scatter3(array_positions(:, 1), array_positions(:, 2), array_positions(:, 3), 100, colors{mod(p-1, length(colors))+1}, markers{mod(p-1, length(markers))+1}, 'filled');

        legend_handles = [legend_handles, h];
        legend_labels{end+1} = placements{p};

        for s = 1:size(source_positions, 1)
            source_position = source_positions(s, :);
            for t_idx = 1:size(all_predictions, 2)
                predicted_position = all_predictions{p, t_idx, s};
                scatter3(predicted_position(1), predicted_position(2), predicted_position(3), 50, colors{mod(p-1, length(colors))+1}, markers{mod(p-1, length(markers))+1});
            end
        end
    end
%}
    h_source = scatter3(source_positions(:, 1), source_positions(:, 2), source_positions(:, 3), 150, 'r', 'd', 'filled');
    legend_handles = [legend_handles, h_source];
    legend_labels{end+1} = 'True Sources';

    legend(legend_handles, legend_labels);
    grid on;
    xlabel('X (cm)'); ylabel('Y (cm)'); zlabel('Z (cm)');
    title('3D Visualization of Microphone Arrangements and Predictions');
    view(3);
    hold off;

    %% Save results
    saveas(gcf, fullfile(output_dir, '3D_visualization.png'));
    saveas(gcf, fullfile(output_dir, '3D_visualization.fig'));

    % Save tables and errors to .mat file
    output_file = fullfile(output_dir, 'simulation_results.mat');
    save(output_file, 'errors', 'T_snr_freq', 'T_freq', 'T_snr', 'T_stats', 'all_predictions');
    disp(['Data saved to ', output_file]);
end
%{
function T_highlighted = highlight_min_with_star(T)
    % Highlight the smallest values in each column
    data = table2array(T); % Extract numeric data
    highlighted_rows = num2cell(data); % Convert data to cell array
    for col = 1:size(data, 2)
        [~, min_idx] = min(data(:, col)); % Find index of the minimum value
        % Append '*' to the smallest value without altering the formatting
        highlighted_rows{min_idx, col} = [num2str(data(min_idx, col)), '*'];
    end
    T_highlighted = cell2table(highlighted_rows, ...
        'VariableNames', T.Properties.VariableNames, ...
        'RowNames', T.Properties.RowNames); % Convert back to table
end

function T_highlighted = highlight_min_with_star(T)
    % Highlight the smallest values in each column
    data = table2array(T); % Extract numeric data
    highlighted_rows = num2cell(data); % Convert data to cell array
    for col = 1:size(data, 2)
        [~, min_idx] = min(data(:, col)); % Find index of the minimum value
        % Append '*' to the smallest value without altering formatting
        highlighted_rows{min_idx, col} = [num2str(data(min_idx, col)), '*'];
    end
    T_highlighted = cell2table(highlighted_rows, ...
        'VariableNames', T.Properties.VariableNames, ...
        'RowNames', T.Properties.RowNames); % Convert back to table
end

function T_highlighted = highlight_min_with_star(T)
    % Highlight the smallest values in each column
    data = table2array(T); % Extract numeric data
    highlighted_rows = data; % Use numeric data array directly
    for col = 1:size(data, 2)
        [~, min_idx] = min(data(:, col)); % Find the minimum value index
        % Append '*' to the smallest value without converting to cell
        highlighted_rows(min_idx, col) = string([num2str(data(min_idx, col)), '*']);
    end

    % Convert highlighted_rows back to table format
    T_highlighted = array2table(highlighted_rows, 'VariableNames', T.Properties.VariableNames, 'RowNames', T.Properties.RowNames);
end

function T_highlighted = highlight_min_with_star(T)
    % Highlight the smallest values in each column
    data = table2array(T); % Extract numeric data
    highlighted_rows = num2cell(data); % Convert to cell array
    for col = 1:size(data, 2)
        [~, min_idx] = min(data(:, col)); % Find the minimum value index
        % Append '*' to the smallest value
        highlighted_rows{min_idx, col} = strcat(num2str(data(min_idx, col)), '*'); % Add star marker while preserving numeric value
    end

    % Convert highlighted_rows back to table format
    T_highlighted = cell2table(highlighted_rows, 'VariableNames', T.Properties.VariableNames, 'RowNames', T.Properties.RowNames);
end

function T_highlighted = highlight_min_with_star(T)
    % Highlight the smallest values in each column
    data = table2array(T); % Extract numeric data
    highlighted_rows = strings(size(data)); % Initialize as strings to avoid brackets
    for col = 1:size(data, 2)
        for row = 1:size(data, 1)
            % Convert each value to string
            highlighted_rows(row, col) = num2str(data(row, col));
        end
        % Find the smallest value and append '*'
        [~, min_idx] = min(data(:, col)); % Find index of the minimum value
        highlighted_rows(min_idx, col) = strcat(highlighted_rows(min_idx, col), '*'); % Append star marker
    end

    % Convert highlighted_rows back to table format
    T_highlighted = array2table(highlighted_rows, 'VariableNames', T.Properties.VariableNames, 'RowNames', T.Properties.RowNames);
end

function T_highlighted = highlight_min_with_star(T)
    % Highlight the smallest values in each column
    data = table2array(T); % Extract numeric data
    highlighted_rows = cell(size(data)); % Initialize as cell array to avoid " "
    for col = 1:size(data, 2)
        for row = 1:size(data, 1)
            % Convert each value to string
            highlighted_rows{row, col} = num2str(data(row, col));
        end
        % Find the smallest value and append '*'
        [~, min_idx] = min(data(:, col)); % Find index of the minimum value
        highlighted_rows{min_idx, col} = strcat(highlighted_rows{min_idx, col}, '*'); % Append star marker
    end

    % Convert highlighted_rows back to table format
    T_highlighted = cell2table(highlighted_rows, 'VariableNames', T.Properties.VariableNames, 'RowNames', T.Properties.RowNames);
end

%}
function T_highlighted = highlight_min_with_star(T)
    % Highlight the smallest values in each column
    data = table2array(T); % Extract numeric data
    highlighted_rows = num2cell(data); % Use cell array for efficient numeric handling

    for col = 1:size(data, 2)
        [~, min_idx] = min(data(:, col)); % Find the smallest value index
        % Append '*' to the smallest value
        highlighted_rows{min_idx, col} = [num2str(data(min_idx, col), '%.5f'), '*'];
    end

    % Convert highlighted_rows back to table format
    T_highlighted = cell2table(highlighted_rows, 'VariableNames', T.Properties.VariableNames, 'RowNames', T.Properties.RowNames);

end