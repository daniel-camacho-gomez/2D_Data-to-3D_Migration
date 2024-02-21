function [D_rat_sim,sigma_sim,desp_sim] = migrationModel(gamma_phi,gamma_theta,v_mean,eta,n_cells,t_assay,Dt_photo,epoch,iter)

for sim=1:n_cells   

    %...Time steps
    Dt_mech    = Dt_photo;     %min
    Dt_mig     = Dt_photo;     %min
    
    t_mig      = Dt_mig;       %counter for migration
    
    t_sim      = 0:Dt_mech:t_assay;   %min         
    
    %-------------------------------------------------------------------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~VARIABLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    %-------------------------------------------------------------------------
    %...cell position           
    x_c =  zeros(length(t_sim),1);     
    y_c =  zeros(length(t_sim),1);      
    z_c =  zeros(length(t_sim),1);     
    
    %...cell velocity
    v_x =  zeros(length(t_sim),1);      
    v_y =  zeros(length(t_sim),1);      
    v_z =  zeros(length(t_sim),1);      
      
    %-------------------------------------------------------------------------
    %~~~~~~~~~~~~~~~~~~~~~~~~Initialization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    %-------------------------------------------------------------------------
    R_c   = 10; %um 
    %...locomotive foce
    F_loc   = 6*v_mean*pi*eta*R_c;
    
    %...Cauchy cdf
    x_0   = 0;
    
    counter_mig = 1; 
    %...initial unit orientation vector
    e_x(counter_mig) = -1+2*rand(1,1);
    e_y(counter_mig) = -1+2*rand(1,1);
    e_z(counter_mig) = 0;
        
    e = [e_x(counter_mig) e_y(counter_mig) e_z(counter_mig)]; 
    e = e/sqrt(e(1)^2+e(2)^2+e(3)^2);
    
    %normalized
    e_x(counter_mig) = e(1);
    e_y(counter_mig) = e(2);
    e_z(counter_mig) = e(3);
    
    %...orientation 
    phi(counter_mig)   = acos(e_z(counter_mig));    
    theta(counter_mig) = atan2(e_y(counter_mig),e_x(counter_mig));
    
    v_x(1:2) =  (1/(6*pi*eta*R_c))*F_loc*e(1);
    v_y(1:2) =  (1/(6*pi*eta*R_c))*F_loc*e(2);
    v_z(1:2) =  (1/(6*pi*eta*R_c))*F_loc*e(3);
    
    %-------------------------------
    %----------Representation------ 
    %-------------------------------      
    [xs,ys,zs] = sphere;     
    xs = xs*R_c;
    ys = ys*R_c;  
    zs = zs*R_c; 
    
    %-------------------------------------------------------------------------
    %~~~~~~~~~~~~~~~~~~~~~~~~~SIMULATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    %-------------------------------------------------------------------------     
    option_solve = optimset('Display','off');
    
    for n=2:length(t_sim)-1
    
        %...Position calculation
        x_c(n+1) = x_c(n) + 0.5*Dt_mech*(3*v_x(n) - v_x(n-1)); 
        y_c(n+1) = y_c(n) + 0.5*Dt_mech*(3*v_y(n) - v_y(n-1)); 
        z_c(n+1) = z_c(n) + 0.5*Dt_mech*(3*v_z(n) - v_z(n-1));
    
        %...direction of migration     
        if t_sim(n) >= t_mig
            %Rotations 
            %...rotation phi
            %inverse transform sampling, range [0.1 0.9]
            sample = 0.1 + 0.8*rand(1,1); 
            %cauchy function
            rot_phi   =  fsolve(@(ang) (1/pi)*atan((ang-x_0)/gamma_phi)+ 1/2 - sample,0,  option_solve); 
            %...rotation theta  
            %inverse transform sampling, range [0.1 0.9]
            sample = 0.1 + 0.8*rand(1,1);   
            %cauchy function
            rot_theta = fsolve(@(ang) (1/pi)*atan((ang-x_0)/gamma_theta)+ 1/2 - sample,0,  option_solve);
            
            %new orientation 
            phi(counter_mig+1)   = phi(counter_mig)   + rot_phi; 
            theta(counter_mig+1) = theta(counter_mig) + rot_theta;
                       
            %new directions
            e_x(counter_mig+1) = sin(phi(counter_mig+1))*cos(theta(counter_mig+1)); 
            e_y(counter_mig+1) = sin(phi(counter_mig+1))*sin(theta(counter_mig+1));     
            e_z(counter_mig+1) = cos(phi(counter_mig+1));
            
            e = [e_x(counter_mig+1) e_y(counter_mig+1) e_z(counter_mig+1)]; 
            e = e/sqrt(e(1)^2+e(2)^2+e(3)^2);
            
            counter_mig = counter_mig + 1;    
            t_mig = t_mig + Dt_mig;
        end   
    
            
        %...Velocity calculation            
        v_x(n+1) =  (1/(6*pi*eta*R_c))*F_loc*e(1);
        v_y(n+1) =  (1/(6*pi*eta*R_c))*F_loc*e(2);
        v_z(n+1) =  (1/(6*pi*eta*R_c))*F_loc*e(3);
    
        
    %------------------------------------
    %----------ANIMATION-----------------
    %------------------------------------
    %             hSurface = surf(xs+x_c(n),ys+y_c(n),zs+z_c(n));
    %                   set(hSurface,'FaceColor',[0 0 1], ...
    %                   'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
    %                   daspect([1 1 1]);
    %             camlight
    % xlim([x_dom(1) x_dom(2)])
    % ylim([y_dom(1) y_dom(2)]) 
    % zlim([z_dom(1) z_dom(2)])
    % 
    % xlabel('$X\;(\mu m)$','Interpreter','Latex','FontSize',15)
    % ylabel('$Y\;(\mu m)$','Interpreter','Latex','FontSize',15)
    % zlabel('$Z\;(\mu m)$','Interpreter','Latex','FontSize',15)
    % pause(0.01)
    % hold off       
    %------------------------------------
    %----------ANIMATION 2D-----------------
    %------------------------------------ 
    %    plot(x_c(n),y_c(n),'ro','Markersize',20)
    %    hold on
    %    quiver(x_c(n),y_c(n),e(1),e(2),'Maxheadsize',20)
    %    ylim([-5 5])
    %    xlim([-5 5])
    %    pause(0.001)
    %    hold off

    end

    %Simulation metrics
    %...directionality ratio for n simulations
    D_cum = 0;
    for n=1:length(t_sim)-1
        d_real = sqrt((x_c(n+1)-x_c(1))^2+(y_c(n+1)-y_c(1))^2+(z_c(n+1)-z_c(1))^2);
        D_cum = D_cum + sqrt((x_c(n+1)-x_c(n))^2 +...
                             (y_c(n+1)-y_c(n))^2 +...
                             (z_c(n+1)-z_c(n))^2);  
        D_rat(n+1,sim) = d_real/D_cum;
         
    end
    
    %...final position
    x_sim(:,sim)   = x_c;
    y_sim(:,sim)   = y_c; 
    z_sim(:,sim)   = z_c; 

end

%...int mean directionalty
D_rat_sim  = trapz((3:length(t_sim))*Dt_mech,mean(D_rat(3:end,:),2));

%...standard deviations
sigma_x    = std(x_sim(end,:));
sigma_y    = std(y_sim(end,:)); 
sigma_sim  = sqrt(sigma_x^2+sigma_y^2); 

%...mean displacements
x_mean = mean(abs(x_sim(end,:))); 
y_mean = mean(abs(y_sim(end,:))); 
desp_sim = sqrt(x_mean^2+y_mean^2);


savedir = fullfile(cd, 'Results');
fname = sprintf('ExampleSim_Epoch_%d_Iter_%d.mat', epoch, iter);
save(fullfile(savedir,fname));
end
        









        
