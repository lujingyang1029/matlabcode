close all;

clear all;

clc;

Iteration = 599; 

T = .1;  

Fx = [1 T; 0 1] ; % state transition matrix

F = blkdiag(Fx, Fx);

Taox = [T^2/2; T]; %input control matrix, this is for the system noise. Considering
%the case only one dimension, Tao is is 2by1 vector. As for the 2-D model,
%Tao should be like as below,

Tao = [Taox; Taox];

% H = [1 0]; % measurement matrix

Nm=100;

%% define main variables

v_mag = .5;  %process noise:

w_mag = [sqrt(3), sqrt(4), sqrt(1)];  %measurement noise pre_set up

R = diag(w_mag.^2); %covariance of measurement noise, in this case, only position 
                %is measuered, so R is a 1by1matrix. 

Qx = v_mag^2 * [T^4/4 T^3/2; T^3/2 T^2]; 

Q = blkdiag(Qx,Qx);

Px=v_mag^2*[1 1/T;1/T 2/(T^2)];

P = blkdiag(Px,Px);

    A=randn(Iteration,1);

    A2=randn(Iteration,1);
    
    w= randn(1,Iteration)*w_mag(1);
    
 
%% initize result variables
% simulate the observation over time
    X_true=[];
        
    Z = []; % real_measurement
    
    X = [430,-3,2438,-15]';
    
%     X=(mvnrnd(X,P,1))';
    
    for t = 1 : Iteration
    
        X= F * X + [v_mag * randn*[(T^2/2); T];v_mag * randn*[(T^2/2); T]];              
            
            X_true=[X_true;X'];
 
    end
    
    


    %Meausrement;

%     target = load('sim1_target_location.mat');
    
%     X_true(:,1)= target.res(:,1);
%     
%     X_true(:,3)= target.res(:,1);
%         
    trajc_x = X_true(:,1);

    trajc_y = X_true(:,3);

%% Simulate the SigIn data at the location (231, 2306)
    dir_vec(:,1) = 231 - X_true(:,1);
    
    dir_vec(:,2) = 2306 - X_true(:,3);
    
    dir_vec(:,3) = sqrt(dir_vec(:,1).^2+ dir_vec(:,2).^2);
    
    dir_vec(:,1)= dir_vec(:,1).*X_true(:,2);
    
    dir_vec(:,2)= dir_vec(:,2).*X_true(:,4);
    
    dir_vec(:,4)= (dir_vec(:,1)+dir_vec(:,2))./dir_vec(:,3);
    
 % The ridial rate is saved in the dir_vec(:,4), [dx,dy]*[vx vy]'/|dx,dy|   
    
    
%% Simulate the DOA estimation;

    theta2_real = atand((trajc_y-2450)./(trajc_x-0));

    theta1_real = atand((trajc_y-1950)./(trajc_x-0));
 
%%    figure;
    



    x_pest = 500./(tand(theta1_real)-tand(theta2_real));

    y_pest = x_pest .* tand(theta1_real) +1950;

    x_pest = x_pest + 0;

%     plot(trajc_x, trajc_y), grid;

%     hold on;

%     plot(x_pest, y_pest,'r*');
    
    Z = [theta1_real theta2_real dir_vec(:,4)];
    
    Z = Z + randn(size(Z))*diag(w_mag);
    
    
    %% Do kalman filtering
    %initize estimation variables
    
    % re-initized state
    
%     X_hat =(mvnrnd(X,P,1))';
    
    X_hat = [430,-3,2438,-15]';
   
    X_pred = [];
    
    predic_var = [];
    
    W_save=[];
    
    Z_est =[];
    
     Z_est_backup=[];
    
    for t = 1:length(X_true)
            
        X_hat = F * X_hat;
            
        P = F * P * F' + Q;
            
        predic_var = [predic_var; P];
        
       % Start of Jacobian matrix for the observation function h(x);
        Deno1 = (X_hat(1)-0)^2+(X_hat(3)-1950)^2;
        
        Nume1x = X_hat(1)-0;
        
        Nume1y = 1950- X_hat(3);
        
        Deno2 = (X_hat(1)-0)^2+(X_hat(3)-2450)^2;
        
        Nume2x = X_hat(1)-0;
        
        Nume2y = 2450- X_hat(3);
        
        % For SigIn the Jacobian row is 
        
        
        
        a = X_hat(1)- 231;
        
        b = X_hat(3) - 2306;
        
        c = a^2 + b^2;
        
        item1 = -X_hat(2)* sqrt(c) + (a*X_hat(2)+b*X_hat(4) )* a/sqrt(c);
        
        item1 = item1/c;
        
        item2 = -a/sqrt(c);
        
        item3 = -X_hat(4)*sqrt(c) + (a*X_hat(2)+b*X_hat(4))*b/sqrt(c);
        
        item3 = item3/c;
        
        item4 = -b/sqrt(c);
        
        H_l= [item1, item2,item3, item4];
              
        H = 180/3.14*[Nume1y/Deno1, 0 , Nume1x/Deno1, 0; Nume2y/Deno2, 0 , Nume2x/Deno2, 0];   
        
        H = [H; H_l];
        
        %End of Jacobian matrix for the observation function h(x);        
          
        W = P*H'/(H*P*H'+R);
        
        W_save=[W_save ; W]; % save the gain matrix
        
        theta2 = atand((X_hat(3)-2450)./(X_hat(1)-0));

        theta1 = atand((X_hat(3)-1950)./(X_hat(1)-0));
        
        s1 = 231 - X_hat(1);
    
        s2 = 2306 - X_hat(3);

        s3 = s1^2+ s2^2;
        
        s4 = s1*X_hat(2) + s2*X_hat(4);
        
        s_est = s4/sqrt(s1^2+s2^2);
        
        Z_est = [theta1; theta2; s_est];
        
        Z_est_backup = [Z_est_backup Z_est];
        
        X_hat = X_hat + W * (Z(t,:)' - Z_est);

        P =  (eye(4)-W*H)*P;
                
        X_pred = [X_pred; X_hat'];

    end
    
    
    figure, plot(theta1_real,'-+'), hold on, plot(theta2_real,'-.s');
    
    plot(Z_est_backup(1,:)','r-<'), plot(Z_est_backup(1,:)', 'k->');
  
    
    figure;
    
    scatter(trajc_x, trajc_y,[],1:length(trajc_x)),hold on;
    
    scatter(X_pred(:,1), X_pred(:,3), [], 1:length(X_pred(:,1))), colorbar;

    save('phase1_2doa_1sigin', 'X_true', 'X_pred','theta1_real','theta2_real','Z_est_backup','dir_vec')

