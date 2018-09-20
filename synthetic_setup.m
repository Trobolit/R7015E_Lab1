RANDOMIZE_AGAIN = 1;

[uDataset1, yDataset1] = textread('Dataset1.txt','%f %f');
[uDataset2, yDataset2] = textread('Dataset2.txt','%f %f');
[uDataset3, yDataset3] = textread('Dataset3.txt','%f %f');

Ad = [1,1;0,1];
Bd = [0;1];
Cd = [1,0.5];
Dd = 0;
Q = [20, 0; 0, 0.1];
x0 = [0;0];

k_max = numel(uDataset1);

x0 = [0;0];

if RANDOMIZE_AGAIN
    w = nan(2, 1000);
    tmp = nan(2,2);
    for k=1:1000
        tmp = normrnd(0,Q);
        w(:,k) = [tmp(1,1); tmp(2,2)];
    end

    v = nan(1000,1);
    for k=1:1000
        v(k) = normrnd(0,4);
    end
    fprintf('\nYou have now resetted the random noise variables data!\n\n');
end
