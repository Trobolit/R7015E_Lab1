%% System Setup
% This system is a pure integrator. The noisless system with given input
% gives a smooth 3rd degree poly looking curve.
Ts = 1;
my_sys = ss(Ad,Bd,Cd,Dd,Ts); %discrete time with sampling period 1s.
my_poles = pole(my_sys);
% Both poles are in 1.

% Choose some reasonable poles
L_poles = [0.4, 0.7]; %[0.7+0.3i; 0.7-0.3i];
Ld = place(Ad', Cd', L_poles)';

% init variables
xhat = nan(2,k_max);
yhat = nan(1,k_max);
ytilde = nan(1,k_max);


%% Luenberger
fprintf('Luenberger:\n');
figure(1);
hold on;
for dataset=1:3

    %initial conditions
    xhat(:,1) = x0;
    yhat(1) = Cd*xhat(:,1);
    switch dataset
        case 1
            yDataset = yDataset1;
            uDataset = uDataset1;
        case 2
            yDataset = yDataset2;
            uDataset = uDataset2;
        case 3
            yDataset = yDataset3;
            uDataset = uDataset3;
    end
    ytilde(1) = yDataset(1) - yhat(1);

    %loop to estimate with Lagrange
    for k=1:k_max-1
        xhat(:,k+1) = Ad*xhat(:,k) + Bd*uDataset(k) + Ld*ytilde(k); % estimate new state
        yhat(k+1) = Cd*xhat(:,k+1);     % calculate yhat for estimated state
        ytilde(k+1) = yDataset(k+1) - yhat(k+1);
    end
    
    %Plot that
    plot(yDataset);
    plot(yhat);
    %print mean square error
    %fprintf('MSE for dataset %d: %d\n',dataset, immse(yhat(2:end),yDataset(1:end-1)') );
    fprintf('MSE for dataset %d: %d\n',dataset, immse(yhat,yDataset') );
     % if you offset these one index performance is 10x better.
     
end

legend('y1','yhat1','y2','yhat2','y3','yhat3');
hold off

%% Stationary Kalman
fprintf('Stationary kalman:\n');
figure(2);
hold on;

% R measurment noise.
% N the way noise propagates through system, design param
% Q Model/state noise
% S cov of noises (Q, R).

% | Q  S |
% | S' R |

% stat kalman init
nu = nan(1,k_max);
%nuhat = nan(2,k_max);
N = eye(2); % Design parameter, the more we trust measurments.
%R1 = N(1); 
%R12 = 0; %assume correlation between w and v is 0.
S = [0;0];
R = 4; %intensity of noise in v. 4 is its standard deviation. ie most values here are below +-4. check: histogram(v)
[P, ~, ~] = dare(Ad', Cd', Q, R);
K = (Ad*P*Cd' + N*S)/( Cd*P*Cd' + R);


for dataset=1:3

    %initial conditions
    xhat(:,1) = x0;
    yhat(1) = Cd*xhat(:,1); 
    nuhat(1) = 0;
    switch dataset
        case 1
            yDataset = yDataset1;
            uDataset = uDataset1;
        case 2
            yDataset = yDataset2;
            uDataset = uDataset2;
        case 3
            yDataset = yDataset3;
            uDataset = uDataset3;
    end
    nu(1) = yDataset(1) - Cd*xhat(:,1);
    

    %loop to estimate with stat kalman
    for k=1:k_max
        % Update current state estimate
        nu(k) = yDataset(k) - Cd*xhat(:,k);
        %nuhat(:,k) = (Q*S/(Cd*P*Cd' + R))*nu(k); % always 0 since S=0.
        xhat(:,k) = xhat(:,k) + K*nu(k);
        yhat(k) = Cd*xhat(:,k);
        if k<k_max
            % Predict next step using updated current state
            xhat(:,k+1) = Ad*xhat(:,k) + Bd*uDataset(k); % + nuhat(:,k);
        end
        
    end
    
    
    %Plot that
    plot(yDataset);
    plot(yhat);
    %print mean square error
    fprintf('MSE for dataset %d: %d\n',dataset, immse(yhat,yDataset') ); %yhat(2:end),yDataset(1:end-1)
     % if you offset these one index performance is 10x better.
 
end

legend('y1','yhat1','y2','yhat2','y3','yhat3');
hold off

%% Non-Stationary Kalman
fprintf('Non-Stationary kalman:\n');
figure(3);
hold on;

% non stat kalman init
N = 1*eye(2); %1 leads to better results than 5 and 0.5
S = [0;0]; %assume correlation between w and v is 0.
R = 4; %intensity of noise in v. 4 is its standard deviation. ie most values here are below +-4. check: histogram(v)
[P0, ~, ~] = dare(Ad', Cd', Q, R);

P = nan(2,2,k_max);
P(:,:,1) = P0;

for dataset=1:3

    %initial conditions
    xhat(:,1) = x0;
    yhat(1) = Cd*xhat(:,1); 
    switch dataset
        case 1
            yDataset = yDataset1;
            uDataset = uDataset1;
        case 2
            yDataset = yDataset2;
            uDataset = uDataset2;
        case 3
            yDataset = yDataset3;
            uDataset = uDataset3;
    end
    

    %loop to estimate with nonstat
    for k=1:k_max
        
        % update current state
        K = P(:,:,k)*Cd'/(Cd*P(:,:,k)*Cd' + R);
        ytilde(k) = yDataset(k) - Cd*xhat(:,k);
        xhat(:,k) = xhat(:,k) + K*ytilde(k); % This is the most accurate state we can have this iteration given now and last.    
        P(:,:,k) = P(:,:,k) - K*Cd*P(:,:,k);
        yhat(k) = Cd*xhat(:,k);
        % Predict
        if k<k_max
            xhat(:,k+1) = Ad*xhat(:,k) + Bd*uDataset(k);
            P(:,:,k+1) = Ad*P(:,:,k)*Ad' + N*Q*N';
        end
        
    end
    
    %Plot that
    plot(yDataset);
    plot(yhat);
    %print mean square error
    fprintf('MSE for dataset %d: %d\n',dataset, immse(yhat,yDataset') ); %immse(yhat(2:end),yDataset(1:end-1)')
 
end

legend('y1','yhat1','y2','yhat2','y3','yhat3');
hold off