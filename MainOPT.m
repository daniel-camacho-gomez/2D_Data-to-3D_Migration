%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Code for predicting 3D migration from 2D data 
% By Daniel Camacho-Gomez,
% Unversity of zaragoza, Spain
% E-mail: dcamacho@unizar.es, danielcamachogomez@hotmail.com
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
clearvars; close all
mkdir Results
%---------------------------------------
%----------EXPERIMENTS------------------ 
%---------------------------------------

%...Experimental assay
t_assay  = 60;  %min duration of the experiment 
Dt_photo = 0.5; %min photo time interval 

%...(uncomment for an example)
%...PBMCs 4 mg/mL
D_rat_exp = 22.32;
sigma_exp = 110.64;
desp_exp  = 81.12;   %um
n_cells   = 75;
v_mean    = 6.80;    %um/min 

% %...PBMCs 6 mg/mL
% D_rat_exp = 19.20;
% sigma_exp = 81.64;
% desp_exp  = 52.40;   %um
% n_cells   = 44;  
% v_mean    = 5.26;    %um/min 

%...CART 4 mg/mL
% D_rat_exp = 20.41;
% sigma_exp = 122.04;
% desp_exp  = 80.17;   %um
% n_cells   = 31;
% v_mean    = 6.53;    %um/min

%...CART 6 mg/mL
% D_rat_exp = 18.11;
% sigma_exp = 70.15;
% desp_exp  = 41.57;   %um
% n_cells   = 49;
% v_mean    = 3.94;    %um/min

%...CXCL12 6 mg/mL
% D_rat_exp = 26.67;
% sigma_exp = 78.42;
% desp_exp  = 51.22;   %um
% n_cells   = 27;
% v_mean    = 4.30;    %um/min 



%-------------------------------
%----------Matrix--------------- 
%-------------------------------        
eta_twofive  = 7.96*1e3*60;   % Pa s --> ug/um min (1e3*60)   2.5 mg/ml
eta_four     = 18.42*1e3*60;  % Pa s --> ug/um min (1e3*60)   4.0 mg/ml
eta_six      = 39.15*1e3*60;  % Pa s --> ug/um min (1e3*60)   6.0 mg/ml

%...select the matrix of the experiment
eta = eta_six;

%-------------------------------
%----------Optimization--------- 
%-------------------------------  
%...steps and iteration data
N_epoch = 500; 
N_iter  = 4;

%....calibration parameters 
gamma_phi_k       = zeros(N_iter,1);
gamma_theta_k     = zeros(N_iter,1);

%...initial seed.
gamma_phi_k(:)    = 0.1;
gamma_theta_k(:)  = 1; 

%...allocation best results
D_rat_best   = 0; 
prev_fitness = 0; 

fitness      = zeros(N_iter,1);
ref          = 1;

D_rat_sim = zeros(N_iter,1); 
sigma_sim = zeros(N_iter,1); 
desp_sim  = zeros(N_iter,1); 

disp("------------------------------------------------------------------------------------")
disp("                                        OPTIMIZING  " );
disp("------------------------------------------------------------------------------------")
tic
for epoch = 1:N_epoch 
    for iter = 1:N_iter
        
         %...perform simulation  
         [D_rat_sim(iter),sigma_sim(iter),desp_sim(iter)] = migrationModel(gamma_phi_k(iter),gamma_theta_k(iter),...
                                            v_mean,eta,n_cells,t_assay,Dt_photo,epoch,iter); 

         %...calculate fitness                                                 
         fitness(iter)   =  1/3*100.^(-(D_rat_sim(iter)-D_rat_exp).^2/(2*(D_rat_exp)^2))+...
                           1/3*100.^(-(sigma_sim(iter)-sigma_exp).^2/(2*(sigma_exp)^2))+...
                           1/3*100.^(-(desp_sim(iter)-desp_exp).^2/(2*(desp_exp)^2));  

    end
    
    [max_fitness(epoch), id_max]  = max(fitness); 
    [total_best, best_epoch]      = max(max_fitness);
    best_pob(epoch) = id_max; 
    
    elap_h =  fix(toc/3600);
    elap_m =  fix(toc/60)-60*fix(toc/3600);
    elap_s =  toc-60*elap_m-3600*elap_h;
    disp("------------------------------------------------------------------------------------")
    disp("------------------------------------------------------------------------------------")
    disp("                                           Epoch:  " + epoch);
    disp("                                      Best Epoch:  " + best_epoch);
    disp("                              Total lowest error:  " + (1-max(max_fitness)));
    disp("                     Previous epoch lowest error:  " + (1-prev_fitness)); 
    disp("                      Current epoch lowest error:  " + (1-max(fitness)));
    disp("                              Total elapsed time:  " + elap_h + ' h ' +elap_m+ ' min ' +elap_s+ ' s');
    disp("------------------------------------------------------------------------------------")
    disp("------------------------------------------------------------------------------------")


    % DIRECT SEARCH METHOD
    prev_fitness  = max(fitness); 
    % select the best two  
    [max1,id1]    = max(fitness); 
    fitness(id1)  = -1000; 
    best_pob(epoch) = id1; 

    [max2,id2] = max(fitness);  
    fitness(id1) = max1; 
  
    %selection of new values   
    gamma_phi_best(epoch)   = gamma_phi_k(id_max);
    gamma_theta_best(epoch) = gamma_theta_k(id_max);
 
    %refine if the error 
    if (1-fitness(id_max)) < 0.05
        ref = 0.5;
    end

    gamma_phi_k(1)    = gamma_phi_best(epoch)   + ref*(1-fitness(id_max))*gamma_phi_best(epoch);
    gamma_theta_k(1)  = gamma_theta_best(epoch) + ref*(1-fitness(id_max))*gamma_theta_best(epoch);
     
    gamma_phi_k(2)    = gamma_phi_best(epoch)   - ref*(1-fitness(id_max))*gamma_phi_best(epoch);
    gamma_theta_k(2)  = gamma_theta_best(epoch) + ref*(1-fitness(id_max))*gamma_theta_best(epoch);
     
    gamma_phi_k(3)    = gamma_phi_best(epoch)   + ref*(1-fitness(id_max))*gamma_phi_best(epoch);
    gamma_theta_k(3)  = gamma_theta_best(epoch) - ref*(1-fitness(id_max))*gamma_theta_best(epoch);
     
    gamma_phi_k(4)    = gamma_phi_best(epoch)   - ref*(1-fitness(id_max))*gamma_phi_best(epoch);
    gamma_theta_k(4)  = gamma_theta_best(epoch) - ref*(1-fitness(id_max))*gamma_theta_best(epoch);
     
end
disp("                                END OF OPTIMIZATION       " );


%-------------------------------------------------------------
% Optimization results
figure
plot(1:length(gamma_phi_best),gamma_phi_best)
xlabel('Epoch','Interpreter','Latex','FontSize',15)
ylabel('$\gamma_{\varphi}$','Interpreter','Latex','FontSize',15)

figure
plot(1:length(gamma_theta_best),gamma_theta_best)
xlabel('Epoch','Interpreter','Latex','FontSize',15)
ylabel('$\gamma_{\theta}$','Interpreter','Latex','FontSize',15)

figure
plot(1:length(max_fitness),(1-max_fitness))
xlabel('Epoch','Interpreter','Latex','FontSize',15)
ylabel('Error','Interpreter','Latex','FontSize',15)

figure
plot(1:length(best_pob),best_pob)
xlabel('Best Epoch','Interpreter','Latex','FontSize',15)
ylabel('Best Iter','Interpreter','Latex','FontSize',15)
%-------------------------------------------------------------


%Save variables
savedir = fullfile(cd, 'Results');
fname   = sprintf('Migration_Example.mat');
save(fullfile(savedir,fname));



 




