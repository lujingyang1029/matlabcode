function [X_pred]= Kalman_filter3(Z,Q,R,F,H,P)
%% Kalman filter for the constant acceleratin lane keeping model

%Z is the measurement including x and y, 
%Q is system process noise
%R is the measurement noise

    X_hat = [Z(1,1),0,0]';
   
    X_pred = [];
    
    predic_var = [];
    
    W_save=[];
    
    Z_est =[];
    
    Z_est_backup=[];
    
    for t = 1:length(Z)
            
        X_hat = F * X_hat;
        
        Z_est = H*X_hat;
            
        P = F * P * F' + Q;
            
        predic_var = [predic_var; P];
        
        W = P*H'/(H*P*H'+R);
        
        W_save=[W_save ; W]; % save the gain matrix
        
        Z_est_backup = [Z_est_backup Z_est];
        
        X_hat = X_hat + W * (Z(t,:)' - Z_est);

        P =  (eye(length(Q))-W*H)*P;
                
        X_pred = [X_pred; X_hat'];

    end