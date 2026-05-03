S0 = 100;
mu = 0.0003;
sigma = 0.01;
T = 252;
num_paths = 1000;

S = zeros(num_paths, T);
S(:,1) = S0;

for i = 1:num_paths
    for t = 2:T
        Z = randn();
        S(i,t) = S(i,t-1) * exp(mu - 0.5*sigma^2 + sigma*Z);
    end
end

figure; hold on;

for i = 1:num_paths
    plot(S(i,:), 'Color', [0 0 1 0.05]);
end

mean_path = mean(S);
plot(mean_path, 'k', 'LineWidth', 2);

final_prices = S(:,end);
p5 = prctile(final_prices, 5);
p95 = prctile(final_prices, 95);

yline(p5, '--r');
yline(p95, '--r');

fill([time(1), time(end), time(end), time(1)], ...
     [p5, p5, p95, p95], ...
     [0.7 0.7 0.7], ...   
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4);

fprintf('5%% percentile: %.5f, 95%% percentile: %.5f]\n', p5, p95);

title('Simulation');