close all;
clear all;
Iteration = 200; 
T = 1;  
F = [1 T; 0 1] ; % state transition matrix
Tao = [T^2/2; T]; %input control matrix
H = [1 0;1 0]; % measurement matrix
Nm=1000;
%% define main variables
X= [1; 1];  %[position; velocity]
X_hat = X;  %x_estimate 
v_mag = .5;  %process noise:
w1_mag = sqrt(3);  %measurement noise:
w2_mag = sqrt(4);
% w1_mag = 3;  %measurement noise:
% w2_mag = 4;
b_mag=sqrt(3000);
R = [w1_mag^2 0; 0 w2_mag^2]; 
Q = v_mag^2 * [T^4/4 T^3/2; T^3/2 T^2]; 
P=v_mag^2*[1 1/T;1/T 2/(T^2)];
c1=w2_mag^2/(w1_mag^2+w2_mag^2);
c2=w1_mag^2/(w1_mag^2+w2_mag^2);
Qk_1=zeros(Iteration,1);
Qk_2=zeros(Iteration,1);
Qk_3=zeros(Iteration,1);
%%
for i=1:Nm %    Monto Carlo
    A2=randn(Iteration,1);
    A3=randn(Iteration,1);
    W1=randn(Iteration,1)*w1_mag;
    W2=randn(Iteration,1)*w2_mag;
    % make three types of combination of injection
    A=randn(Iteration,1); %produce the bias information
    for j=1:2
        Qk=zeros(Iteration,1);
        MC=0;        
        if j==1
            b1=A*b_mag*c1/(sqrt(c1^2+c2^2));%optimal solution
            b2=A*b_mag*c2/(sqrt(c1^2+c2^2));
        elseif j==2
            b1=A3*b_mag*sqrt(2)/2;
            b2=A2*b_mag*sqrt(2)/2;       
        elseif j==3 %put all to the compromised sensor whose std is smaller  
            b1=A*b_mag;
            b2=A*b_mag*0;      
        end 

% initize result variables
% simulate the observation over time
    X_true=[];
    X1_val = []; % real_locaion
    X2_val = []; % real_velocity
    Z = []; % real_measurement
    X=(mvnrnd(X,P,1))'
    for t = 1 : T: Iteration
        X= F * X + v_mag * A2(t)*[(T^2/2); T]; 
% Generate the measurement       
        if t==Iteration/2
%           y= H* X + w_mag * randn +b_mag * randn;
            y = H * X + [W1(t); W2(t)]+[c1*b1(t); c2*b2(t)];
        else
            y = H * X + [ W1(t); W2(t)];

        end
        X1_val = [X1_val; X(1)];
        X2_val = [X2_val; X(2)];
        X_true=[X_true;X'];
        Z = [Z y];   
    end
    %% Do kalman filtering
    %initize estimation variables
    X2_estimate = []; % Velocity estimate
    X= [0; 0]; % re-initized state
    X_hat = X;
    X_predict_state = [];
    X1_estimate = []; %  Position estimate
    predic_var = [];
    W_save=[];
    Sum=[];
    for t = 1:length(X1_val)
            X_hat = F * X_hat;
%           X_predict_state = [X_predict_state; X_hat(1)] ;% Save the information about the location
            %predict next covariance
            P = F * P * F' + Q;
            predic_var = [predic_var; P] ;
            % predicted measurement covariance
            % Kalman Gain
%             eee=H*P*H'+R;
            W = P*H'*inv(H*P*H'+R);
            W_save=[W_save ; W]; % save the gain matrix
            % Update the state estimate.
            X_hat = X_hat + W * (Z(:,t) - H * X_hat);
            % update covariance estimation.
%             Sum=Sum+(X_hat-[X_val(t);X_dot(t)])'*inv(X)*(X_hat-[X_val(t);X_dot(t)]);
%             Qk=[Qk;Sum];
%             Sum=0;
            P =  (eye(2)-W*H)*P;
            Sum=[Sum;(X_hat-(X_true(t,:))')'*inv(P)*(X_hat-(X_true(t,:))')];
        %Store for plotting      
        X1_estimate = [X1_estimate; X_hat(1)];
        X2_estimate = [X2_estimate; X_hat(2)];
    end
    Qk=Sum+Qk;
    if j==1
        Qk_1=Qk_1+Qk;
    elseif j==2
        Qk_2=Qk_2+Qk;
    elseif j==3
        Qk_3=Qk_3+Qk;
    end
 end
        
end

    figure(1);
    h = gca; % 
    set(h,'FontSize',14); % 
    tt = 1 : T : Iteration;
    tt=95:T:110;


        Qk=Qk_1(95:110);
        plot(tt,Qk,'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',8)
        hold on;
               Qk=Qk_2(95:110);
        plot(tt,Qk,'--go','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','y',...
            'MarkerSize',8)

%         Qk=Qk_3(90:120);
%          plot(tt,Qk,'--bd','LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor','b',...
%             'MarkerSize',8)

    xlabel('Iteration Number k','FontSize',14),ylabel('q_k');
    fprintf(1,'Hope it will work');
    grid;

    legend('Optimal Solution','Equally distributed','All on sensor1')
