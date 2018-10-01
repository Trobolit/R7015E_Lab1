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


title('Luenberger');
legend('x_w','theta','xwhat','thetahat');
hold off
figure(11);
hold on;
plot(ytilde(1,:));
plot(ytilde(2,:));
legend('xw pred err', 'theta pred err');
title('Luenberger pred err');
hold off;

%% Stationary Kalman
free_code;
fprintf('Stationary kalman:\n');

Q = eye(4);

% stat kalman init
nu = nan(2,k_max);
%nuhat = nan(2,k_max);
N = [0.3, 0,0,0;
     0, 0.1, 0,0;
     0, 0, 0.003, 0;
     0, 0, 0, 0.8];

S = zeros(4,2); %assume correlation between w and v is 0.
R = zeros(2,2);
R(1,1) = 1*cov(yDataset(1,:));
R(2,2) = 0.4*cov(yDataset(2,:));


[P, ~, ~] = dare(Ad', Cd', N*Q*N', R);
K = (Ad*P*Cd' + N*S)/( Cd*P*Cd' + R);

%initial conditions
xhat(:,1) = x0;
yhat(:,1) = Cd*xhat(:,1); 
nuhat(1) = 0;
nu(:,1) = yDataset(:,1) - Cd*xhat(:,1);

estimation_error=nan(2,k_max);
prediction_error=nan(2,k_max);


%loop to estimate with stat kalman
for k=1:k_max
    % Update current state estimate
    nu(:,k) = yDataset(:,k) - Cd*xhat(:,k);
    %nuhat(:,k) = (Q*S/(Cd*P*Cd' + R))*nu(k); % always 0 since S=0.
    xhat(:,k) = xhat(:,k) + K*nu(:,k);
    yhat(:,k) = Cd*xhat(:,k);
    estimation_error(:,k) = yDataset(:,k) - yhat(:,k);
    if k<k_max
        % Predict next step using updated current state
        xhat(:,k+1) = Ad*xhat(:,k) + Bd*uDataset(k); % + nuhat(:,k);
        prediction_error(:,k+1) = yDataset(:,k) - Cd*xhat(:,k+1);
    end

end

%clf(2);
figure(2);
hold on;
%Plot that
plot(yDataset(1,:));
plot(yDataset(2,:));
plot(yhat(1,:));
plot(yhat(2,:));
%print mean square error
%fprintf('MSE for x_w %d: %d\n',dataset, immse(yhat(1,:),yDataset(1,:)) );
%fprintf('MSE for theta %d: %d\n',dataset, immse(yhat(2,:),yDataset(2,:)) );
 % if you offset these one index performance is 10x better.

title('stationary');
legend('x_w','theta','xwhat','thetahat');
%legend('x_w','xwhat');
%axis([6150,6300,0.562,0.568]);
hold off

%{
clf(222)
figure(222)
hold on;
%plot(yDataset(1,:));
plot(yDataset(2,:));
%plot(yhat(1,:));
plot(yhat(2,:));
legend('theta','thetahat');
%axis([3240,3380,-0.055,-0.048]);
hold off;
%}


figure(22);
hold on;
plot(prediction_error(1,:));
plot(prediction_error(2,:));
plot(estimation_error(1,:));
plot(estimation_error(2,:));
title('Stationary Kalman errors');
legend('x, pred err', 'theta, pred error', 'x, est err', 'theta, est err');
hold off;

%print mean square error
fprintf('estimation MSE for x_w %d: %d\n',dataset, immse(estimation_error(1,:),zeros(1,k_max)) );
fprintf('estimation MSE for theta %d: %d\n',dataset, immse(estimation_error(2,:),zeros(1,k_max)) );
%print mean square error
fprintf('prediction MSE for x_w %d: %d\n',dataset, immse(prediction_error(1,2:end),zeros(1,k_max-1)) );
fprintf('prediction MSE for theta %d: %d\n',dataset, immse(prediction_error(2,2:end),zeros(1,k_max-1)) );


%% Non-Stationary Kalman
fprintf('Non-Stationary kalman:\n');
figure(3);
hold on;

estimation_error=nan(2,k_max);
prediction_error=nan(2,k_max);

% non stat kalman init
%N = eye(4); %1 leads to better results than 5 and 0.5
N = [0.3, 0,0,0;
     0, 0.1, 0,0;
     0, 0, 0.003, 0;
     0, 0, 0, 0.8];
R(1,1) = 1*cov(yDataset(1,:));
R(2,2) = 0.4*cov(yDataset(2,:));

[P0, ~, ~] = dare(Ad', Cd', Q, R);

P = nan(4,4,k_max);
P(:,:,1) = P0;

%initial conditions
xhat(:,1) = x0;
yhat(:,1) = Cd*xhat(:,1); 

%loop to estimate with nonstat
for k=1:k_max

    % update current state
    K = P(:,:,k)*Cd'/(Cd*P(:,:,k)*Cd' + R);
    ytilde(:,k) = yDataset(:,k) - Cd*xhat(:,k);
    xhat(:,k) = xhat(:,k) + K*ytilde(:,k); % This is the most accurate state we can have this iteration given now and last.    
    P(:,:,k) = P(:,:,k) - K*Cd*P(:,:,k);
    yhat(:,k) = Cd*xhat(:,k);
    estimation_error(:,k) = yDataset(:,k) - yhat(:,k);
    % Predict
    if k<k_max
        xhat(:,k+1) = Ad*xhat(:,k) + Bd*uDataset(k);
        P(:,:,k+1) = Ad*P(:,:,k)*Ad' + N*Q*N';
        prediction_error(:,k+1) = yDataset(:,k) - Cd*xhat(:,k+1);
    end

end

%plot that
plot(yDataset(1,:));
plot(yDataset(2,:));
plot(yhat(1,:));
plot(yhat(2,:));

title('non stationary');
legend('x_w','theta','xwhat','thetahat');
hold off

figure(33);
hold on;
plot(prediction_error(1,:));
plot(prediction_error(2,:));
plot(estimation_error(1,:));
plot(estimation_error(2,:));
title('non-Stationary Kalman errors');
legend('x, pred err', 'theta, pred error', 'x, est err', 'theta, est err');
hold off;

%print mean square error
fprintf('estimation MSE for x_w %d: %d\n',dataset, immse(estimation_error(1,:),zeros(1,k_max)) );
fprintf('estimation MSE for theta %d: %d\n',dataset, immse(estimation_error(2,:),zeros(1,k_max)) );
%print mean square error
fprintf('prediction MSE for x_w %d: %d\n',dataset, immse(prediction_error(1,2:end),zeros(1,k_max-1)) );
fprintf('prediction MSE for theta %d: %d\n',dataset, immse(prediction_error(2,2:end),zeros(1,k_max-1)) );


%% close stuff
close(1);
close(11);
