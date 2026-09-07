clear;
close all;

phi1Data = load('./phi1.m');
phi3Data = load('./phi3.m');

gridSize = size(phi1Data, 2);
snapshotIndices = [3, 21, 101];
panelPositions = [
    0.141904761904762, 0.6314, 0.722023736862908, 0.293512950793925;
    0.141904761904762, 0.3123, 0.722023736862908, 0.293512950793925;
    0.141904761904762, 0.0068, 0.722023736862908, 0.293512950793925
];

r = linspace(0, 1, gridSize);
z = linspace(0, 2, 2 * gridSize);
[R, Z] = meshgrid(r, z);

solidField = phi3Data(1:gridSize, :);
solidField = [flipud(solidField); solidField];

fig = figure(1);
set(fig, 'Position', [1169, 817, 391, 507]);
colormap gray;

for panel = 1:numel(snapshotIndices)
    rows = (snapshotIndices(panel) - 1) * gridSize + (1:gridSize);
    phaseSnapshot = phi1Data(rows, :);
    phaseField = [flipud(phaseSnapshot); phaseSnapshot];

    subplot('Position', panelPositions(panel, :));
    contourf(R, Z, phaseField, [0.5, 0.5], ...
        'FaceColor', 'y', 'EdgeColor', 'k', 'LineWidth', 1, 'LineStyle', '-');
    hold on;
    contourf(R, Z, solidField, [0.5, 0.5], ...
        'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'k', 'LineWidth', 1, 'LineStyle', '-');
    axis image;
    axis([0, 1, 0, 2]);
    xticks([]);
    yticks([]);
    view(-90, 90);
end

print(fig, 'fig2c.pdf', '-dpdf');
