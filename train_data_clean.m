data = load('train_data.mat');
NN_samples = data.NN_samples;
NN_targets = data.NN_targets;

des_samples = 50;

% Get the number of samples for each letter
u = blanks(4*26);
d = u;
ns = 3; % spaces for each char display
for c = 1:26
    samples = find(NN_targets(c,:) == 1);
    n_samples = length(samples);
    if n_samples>des_samples
        % Take some samples out
        n_out = n_samples-des_samples;
        out = randsample(samples,n_out);
        NN_targets(:,out) = [];
        NN_samples(:,out) = [];
    end
    
    u(1+ns*(c-1):ns*c) = sprintf([blanks(ns-1),'%c'],char(96+c));
    d(1+ns*(c-1):ns*c) = sprintf(['%',num2str(ns),'d'], n_samples);
end
% Show samples per letter
disp(u); disp(d);

for c = 1:26
    % Show training samples per letter
    samples = find(NN_targets(c,:) == 1);
    n_samples = length(samples);
    gridn = 7;
    nf = ceil(n_samples/gridn^2);
    for f = 1:nf
        h = figure();
        for i = 1:min(gridn^2,n_samples-(f-1)*gridn^2)
            l = samples(i+(f-1)*gridn^2);
            subplot(gridn,gridn,i)
            imshow(reshape(NN_samples(:,l), 11, 11));
        end
%         uiwait(h) % Prevent all other processes from starting until closed.
    end
end

% Save data
save('train_data.mat','NN_samples','NN_targets');