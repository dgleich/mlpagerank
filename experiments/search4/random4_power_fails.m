%% init with
% found = {};
% save 'random-34.mat' found

rng('shuffle')
while 1
    R = randi([1,2],3,27)-1;
    [Rs,niter] = random_trial(R,0.99);
    if numel(Rs) > 0
        printf('Found example!\n');
        load 'random-34.mat'
        found = [found; R];
        save 'random-34.mat' found
    end
end

