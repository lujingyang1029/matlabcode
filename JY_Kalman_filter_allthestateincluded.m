%% For position sensor only case. 
%  This is for the dependent case. 

close all;

clear all;

clc;

Iteration = 200; 

T = .5;  

% Fx = [1 T; 0 1] ; % state transition matrix

Fx = [1 T 1/2*T^2; 0 1 T; 0 0 1];

F = blkdiag(Fx, Fx);

Taox = [T^2/2; T ;1]; %input control matrix, this is for the system noise. Considering
%the case only one dimension, Tao is is 2by1 vector. As for the 2-D model,
%Tao should be like as below,

Taoy = [T^2/2; T ;1];

H = [1 0 0 0 0 0 ;0 0 0 1 0 0 ]; % measurement matrix

%% define main variables

v_mag1 = 0.01;  %System process noise:

v_mag2 = 0.01; %System process noise;
              
Qx1 = v_mag1^2 * [T^4/4 T^3/2 T^2/2; T^3/2 T^2 T; T^2/2 T 1]; 

Qx2 = v_mag2^2 * [T^4/4 T^3/2 T^2/2; T^3/2 T^2 T; T^2/2 T 1]; 

Q = blkdiag(Qx1,Qx2)

Px1=v_mag1^2*[1 1/T 1/(T^2); 1/T 2/(T^2) 3/(T^3);1/(T^2) 3/(T^3) 6/(T^4) ]

Px2=v_mag2^2*[1 1/T 1/(T^2); 1/T 2/(T^2) 3/(T^3);1/(T^2) 3/(T^3) 6/(T^4) ]

P = blkdiag(Px1,Px2)

w_mag = [sqrt(3), sqrt(4)];  %measurement noise pre_set up

R = diag(w_mag.^2) %covariance of measurement noise, in this case, only position 
                %is measuered, so R is a 1by1matrix. 

 %% initize result variables
% simulate the observation over time
    X_true=[];
        
    Z = []; % real_measurement
    
    X = [100,20,4, 100,30,-1]';
    
%     X=(mvnrnd(X,P,1))';
    
    for t = 1 : Iteration
    
        X= F * X + [v_mag1 * randn*[(T^2/2); T; 1];v_mag2 * randn*[(T^2/2); T; 1]];              
            
            X_true=[X_true;X'];
    end
    
    figure, plot(X_true(:,1),X_true(:,4));
    
    Z = X_true * H';
    
    Z = Z + randn(size(Z))*diag(w_mag);
    
     load JY_xy.mat
    %% Do kalman filtering
    %initize estimation variables
    
    % re-initized state
    
%   X_hat =(mvnrnd(X,P,1))';
    
    X_hat = [Z(1,1),0,0,Z(1,2),0,0]';
   
    X_pred = [];
    
    predic_var = [];
    
    W_save=[];
    
    Z_est =[];
    
    Z_est_backup=[];
    
    for t = 1:length(X_true)
            
        X_hat = F * X_hat;
        
        Z_est = H*X_hat;
            
        P = F * P * F' + Q;
            
        predic_var = [predic_var; P];
        
        W = P*H'/(H*P*H'+R);
        
        W_save=[W_save ; W]; % save the gain matrix
        
        Z_est_backup = [Z_est_backup Z_est];
        
        X_hat = X_hat + W * (Z(t,:)' - Z_est);

        P =  (eye(6)-W*H)*P;
                
        X_pred = [X_pred; X_hat'];

    end
    

    figure;
    
%     scatter(X_true(:,1), X_true(:,4),[],1:length(X_true)),hold on;
    plot(X_true(:,1), X_true(:,4)),hold on;
    
    scatter(X_pred(:,1), X_pred(:,4),[],1:length(X_pred)), colorbar;
    title('Target Tracking')
    xlabel('x'),ylabel('y');
    legend('Ground Truth','Tracking Estimation');
    a = gca;
    a.FontSize = 15;
    grid;
    a.LineWidth= 1.5
    
    figure;
    
    plot(X_true(:,1),'-go'), hold on;
    
    plot(X_pred(:,1),'-.rh')

    title('x');
    
    figure;
    
    plot(X_true(:,4),'-go'), hold on;
    
    plot(X_pred(:,4),'-.rh')

    title('y');
    
    save('JY_xy','X_true','Z');

