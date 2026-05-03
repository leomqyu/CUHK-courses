num_trials = 10000;
win_count = 0;
rounds_total = 0;

for i = 1:num_trials
    balance = 50;
    rounds = 0;
    
    while balance > 0 && balance < 100
        if rand() < 0.48
            balance = balance + 1;
        else
            balance = balance - 1;
        end
        rounds = rounds + 1;
    end
    
    if balance == 100
        win_count = win_count + 1;
    end
    
    rounds_total = rounds_total + rounds;
end

prob_win = win_count / num_trials;
avg_rounds = rounds_total / num_trials;

fprintf('(a) Probability of reaching 100: %.4f\n', prob_win);
fprintf('(b) Average rounds before stopping: %.2f\n', avg_rounds);