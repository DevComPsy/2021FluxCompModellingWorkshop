
function [blocks, choices, rewards] = sim_null(task_struct)

%get task structure data
num_blocks = task_struct.num_blocks;
num_block_trials = task_struct.num_block_trials;
reward_probs = task_struct.reward_probs;


%initialize choices
choices = zeros(num_blocks * num_block_trials, 1);
rewards = zeros(num_blocks * num_block_trials, 1);
blocks = zeros(num_blocks * num_block_trials, 1);

%loop through blocks and trials
for block = 1:num_blocks
    choice_probs = [.5 .5]; % Initialize values for each subject at beginning of block
    
    %loop through trials
    for trial = 1:num_block_trials
        
        % Coin flip to determine which bandit to choose
        if rand(1) < choice_probs(1)
            choice = 1;
        else
            choice = 2;
        end
        
        %Coin flip to determine reward outcome
        if rand(1) < reward_probs(choice)
            reward = 1;
        else
            reward = 0;
        end
    
    %save trial information
    choices(trial + (block - 1) * num_block_trials) = choice;
    rewards(trial + (block - 1) * num_block_trials) = reward;
    blocks(trial + (block - 1) * num_block_trials) = block;
    end
end

end





