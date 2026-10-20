clear;
array = [1 0 0 0 0;
         0 0 1 0 0;
         0 0 0 1 0;
         0 1 0 0 0;
         0 0 0 0 1];
% Find positions of active elements
[row, col] = find(array == 1);

%% 1. 2D Array Geometry Visualization with Grid
figure;

subplot(2,2,1);  % First subplot (top-left)

% Add 5x5 grid overlay first
hold on;
for i = 0.5:1:5.5
    plot([0.5, 5.5], [i, i], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);  % Light gray horizontal lines
    plot([i, i], [0.5, 5.5], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);  % Light gray vertical lines
end

% Plot the dots on top of the grid
scatter(col, row, 600, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'LineWidth', 2);

% Aesthetic settings
title('2D Array Geometry');
xlabel('X-axis (Columns)');
ylabel('Y-axis (Rows)');
axis([0.5 5.5 0.5 5.5]);
axis square;
grid on;

% Hide tick labels
set(gca, 'XTick', [], 'YTick', []);
box on;  % Add border around plot

%% 2. Autocorrelation Plot
% Autocorrelation calculation
autocorr_result = xcorr2(array);

% Prepare meshgrid for plotting
[x, y] = meshgrid(-(size(autocorr_result, 2)/2):(size(autocorr_result, 2)/2 - 1), ...
                  -(size(autocorr_result, 1)/2):(size(autocorr_result, 1)/2 - 1));

subplot(2,2,2);  % Second subplot (top-right)
caxis([0, 3]);
colormap(turbo);  % Hot colormap for better red intensity
surf(x, y, autocorr_result, 'EdgeColor', 'none'); % Smooth surface plot

% Enhanced colormap for strong red contrast
colorbar;
shading interp;
lighting phong;
camlight left;
title('Array Autocorrelation');
xlabel('i');
ylabel('j');
zlabel('\chi_c(i,j)');
view(45, 30);  % Adjust view for 3D effect
axis tight;

%% 3. Distinct Displacement Vectors Visualization with Grid
subplot(2,2,3);  % Third subplot (bottom-left)

% Add 5x5 grid overlay first
hold on;
for i = 0.5:1:5.5
    plot([0.5, 5.5], [i, i], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);  % Light gray horizontal lines
    plot([i, i], [0.5, 5.5], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);  % Light gray vertical lines
end

% Draw displacement vectors between all pairs of elements
for i = 1:length(row)
    for j = 1:length(row)
        % Draw arrow between elements
        quiver(col(i), row(i), col(j)-col(i), row(j)-row(i), 0, ...
               'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', 'r');
    end
end

% Plot the dots on top of arrows and grid
scatter(col, row, 600, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'LineWidth', 2);

% Aesthetic settings
title('Distinct Displacement Vectors');
xlabel('Column Index');
ylabel('Row Index');
axis([0.5 5.5 0.5 5.5]);
axis square;
grid on;

% Hide tick labels
set(gca, 'XTick', [], 'YTick', []);
box on;

%% 4. Enhanced Array Element Density Map (Heatmap with Strong Red)
subplot(2,2,4);  % Fourth subplot (bottom-right)
imagesc(array);  % Display the array as an image

% Enhanced colormap with strong red for intensity
colormap(parula);
caxis([0, max(array(:)) * 1.5]);  % Stretch color limits to intensify red
colorbar;

% Aesthetic settings
title('Enhanced Element Density Map');
xlabel('Columns');
ylabel('Rows');
axis equal;
set(gca, 'YTick', [], 'XTick', []);  % Hide ticks for clean look
box on;

%% Adjust the figure layout
set(gcf, 'Position', [100, 100, 1000, 1000]);  % Resize the figure window for square layout

%% Microphone Placement Configurations
scaling_factor = 3;
x_offset = -10;
y_offset = -4;

costas_array = [
     1, 1;
     2, 4;
     3, 2;
     4, 3;
     5, 5
] * scaling_factor + [x_offset, y_offset];

costas_array_3D = [[costas_array, zeros(size(costas_array, 1), 1)]; [-12, 2.3, 10]];
    markers = {'o', '^', 's', 'd', 'p', 'h', 'x', '+', 'v','<',};
%    colors = {'b', 'm', 'g', 'c', 'y', 'k', 'r', 'b', 'm', 'g'};
colors = {'b', 'm', 'g', 'c', 'y', 'k', 'r', ...  % Standard colors
          [0.5, 0.2, 0.8], [0.9, 0.3, 0.3], [0.3, 0.7, 0.2], ... % Custom RGB colors
          [0.1, 0.5, 0.9], [0.6, 0.2, 0.5], [0.9, 0.7, 0.1], ...
          [0.3, 0.3, 0.3], [0.7, 0.4, 0.1], [0.2, 0.8, 0.8]};
legend_handles = [];
legend_labels = {};
figure;
subplot(1,2,1);
h = scatter3(costas_array(:, 1), costas_array(:, 2), 0, 200, ...
                 'b', 'o','MarkerFaceColor', 'b','MarkerEdgeColor', 'k','LineWidth', 0.5); % Explicitly define edge color
subplot(1,2,2);
    h = scatter3(costas_array_3D(:, 1), costas_array_3D(:, 2), costas_array_3D(:, 3), 200, ...
                 'b', 'o','MarkerFaceColor', 'b','MarkerEdgeColor', 'k','LineWidth', 0.5); % Explicitly define edge color
                 
    legend_handles = [legend_handles, h];
