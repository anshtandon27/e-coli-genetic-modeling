%Import fluorescence data
filename = 'Data.xlsx';
sheetname = 'Sheet2';
data = table2array(readtable(filename, 'Sheet', sheetname));

conc = data(:, 1)*10e6;
avg = data(:, 8:10);
stdev = data(:, 5:7);
AHL = logspace(-3, 5, 100);

%Initialize parameters of genetic circuit
paramInit = [0.1   5.0000e-01   2.3100e-02   1.3000e-05   5.0000e-02   2.0000e-01   2.0000e+00   4.0000e-04];

%Iterate through each strain data
for i = 1:3
    params = optimizer(conc, avg(:, i), paramInit); %Optimizer code
    disp(params);
    modelOutput = SyntheticBio(conc, params);
    RMSE = sqrt(mean((modelOutput-avg(:, i)).^2));
    optimOutput = SyntheticBio(AHL, params);
    figure;
    ax = gca;
    hold on
    errorbar(conc, avg(:, i), stdev(:, i), 'o', 'MarkerSize', 5, 'CapSize', 10);
    plot(AHL, optimOutput);
    set(gca, 'Xscale', 'log');
    xlabel('Concentration (uM)');
    ylabel('Normalized GFP Concentration')
    title(['Strain ' num2str(i)])
    ylim([-0.2, 1.2])
    legend('Experimental', 'Model', 'Location', 'southeast')
    text(max(conc)*0.005, max(avg(:, i))*0.05, ['RMSE: ' num2str(RMSE, '%.3f')], 'FontSize', 10, 'FontWeight', 'bold');
    grid on;
    axis square;
    exportgraphics(ax, sprintf('%d.jpg', i),"Resolution",300)
end




