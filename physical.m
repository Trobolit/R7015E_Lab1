load Dataset4.mat
free_code;
C = [1,0,0,0;0,0,1,0];
my_sys = c2d(ss(A,B,C,D),fSamplingPeriod);

Ad = my_sys.A;
Bd = my_sys.B;
Cd = my_sys.C;
Dd = my_sys.D;

my_poles = pole(my_sys);

Kd; % remainder to me that we have this.

% Choose some reasonable poles
L_poles = exp(Kd*fSamplingPeriod*4);
Ld = place(Ad', Cd', L_poles)';

k_max = numel(afTimes);
x0 = [0;0;0;0];

% init variables
xhat = nan(4,k_max);
yhat = nan(2,k_max);
ytilde = nan(2,k_max);


%% Luenberger
fprintf('Luenberger:\n');
figure(1);
hold on;

%initial conditions
xhat(:,1) = x0;
yhat(:,1) = Cd*xhat(:,1);

yDataset = [aafProcessedInformation(MEASURED_X_W_INDEX,:);aafProcessedInformation(MEASURED_THETA_B_INDEX,:)];
uDataset = aafProcessedInformation(U_INDEX,:);

ytilde(:,1) = yDataset(:,1) - yhat(:,1);

%loop to estimate with Luenberger
for k=1:k_max-1
    xhat(:,k+1) = Ad*xhat(:,k) + Bd*uDataset(k) + Ld*ytilde(:,k); % estimate new state
    yhat(:,k+1) = Cd*xhat(:,k+1);     % calculate yhat for estimated state
    ytilde(:,k+1) = yDataset(:,k+1) - yhat(:,k+1);
end

%Plot that
plot(yDataset(1,:));
plot(yDataset(2,:));
plot(yhat(1,:));
plot(yhat(2,:));
%print mean square error
fprintf('MSE for x_w %d: %d\n',dataset, immse(yhat(1,:),yDataset(1,:)) );
fprintf('MSE for theta %d: %d\n',dataset, immse(yhat(2,:),yDataset(2,:)) );
 % if you offset these one index performance is 10x better.



legend('x_w','theta','xwhat','thetahat');
hold off

%% Stationary Kalman
fprintf('Stationary kalman:\n');
figure(2);
hold on;

Q = eye(4);

% stat kalman init
nu = nan(2,k_max);
%nuhat = nan(2,k_max);
N = eye(4); % why 5? Design parameter? The higher the more we tryst y?
R1 = cov(yDataset(1,:));
R12 = 0; %assume correlation between w and v is 0.
R2 = cov(yDataset(2,:)); %intensity of noise in v. 4 is its standard deviation. ie most values here are below +-4. check: histogram(v)
R = [R1,R12;R12,R2];
[P, ~, ~] = dare(Ad', Cd', N*Q*N', R);
Ktilde = P*Cd'/(Cd*P*Cd' + R2);

%initial conditions
xhat(:,1) = x0;
yhat(:,1) = Cd*xhat(:,1); 
nuhat(1) = 0;
nu(:,1) = yDataset(:,1) - Cd*xhat(:,1);


%loop to estimate with stat kalman
for k=1:k_max
    % Update current state estimate
    nu(:,k) = yDataset(:,k) - Cd*xhat(:,k);
    %nuhat(:,k) = (R1*R12/(Cd*P*Cd' + R2))*nu(k); % always 0 since R12=0.
    xhat(:,k) = xhat(:,k) + Ktilde*nu(:,k);
    yhat(:,k) = Cd*xhat(:,k);
    if k<k_max
        % Predict next step using updated current state
        xhat(:,k+1) = Ad*xhat(:,k) + Bd*uDataset(k); % + nuhat(:,k);
    end

end


%Plot that
plot(yDataset(1,:));
plot(yDataset(2,:));
plot(yhat(1,:));
plot(yhat(2,:));
%print mean square error
fprintf('MSE for x_w %d: %d\n',dataset, immse(yhat(1,:),yDataset(1,:)) );
fprintf('MSE for theta %d: %d\n',dataset, immse(yhat(2,:),yDataset(2,:)) );
 % if you offset these one index performance is 10x better.

legend('x_w','theta','xwhat','thetahat');
hold off

%% Non-Stationary Kalman
fprintf('Non-Stationary kalman:\n');
figure(3);
hold on;

% non stat kalman init
N = eye(4); %1 leads to better results than 5 and 0.5

[P0, ~, ~] = dare(Ad', Cd', N*Q*N', R);

P = nan(4,4,k_max);
P(:,:,1) = P0;

%initial conditions
xhat(:,1) = x0;
yhat(:,1) = Cd*xhat(:,1); 

%loop to estimate with nonstat
for k=1:k_max

    % update current state
    K = P(:,:,k)*Cd'/(Cd*P(:,:,k)*Cd' + R2);
    ytilde(:,k) = yDataset(:,k) - Cd*xhat(:,k);
    xhat(:,k) = xhat(:,k) + K*ytilde(:,k); % This is the most accurate state we can have this iteration given now and last.    
    P(:,:,k) = P(:,:,k) - K*Cd*P(:,:,k);
    yhat(:,k) = Cd*xhat(:,k);
    % Predict
    if k<k_max
        xhat(:,k+1) = Ad*xhat(:,k) + Bd*uDataset(k);
        P(:,:,k+1) = Ad*P(:,:,k)*Ad' + N*Q*N';
    end

end

%plot that
plot(yDataset(1,:));
plot(yDataset(2,:));
plot(yhat(1,:));
plot(yhat(2,:));
%print mean square error
fprintf('MSE for x_w %d: %d\n',dataset, immse(yhat(1,:),yDataset(1,:)) );
fprintf('MSE for theta %d: %d\n',dataset, immse(yhat(2,:),yDataset(2,:)) );

legend('x_w','theta','xwhat','thetahat');
hold off