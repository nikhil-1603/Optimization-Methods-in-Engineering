clear all ;
clc;

%Q5 Zakharov  function
function fun_val = Objective_Fun(x,N) 
    sum1 = 0 ;
    for i=1:N
        sum1 = sum1 + (x(i))^2 ;
    end
    sum2 = 0 ;
    for i=1:N
        sum2 = sum2 + 0.5*i*x(i) ;
    end 
    sum3 = 0 ;
    for i=1:N
        sum3 = sum3 + 0.5*i*x(i) ;
    end     
    fun_val = sum1 + (sum2)^2  + (sum3)^4 ;
end 


% step 1 ##########################################################
N = 2;
A = -5;
B = 10 ;
eps = 0.001 ;
for i=1:N 
   x_0(i) = A + (B-A)*rand  ; 
end 
%x_0 = [-4.0279,4.7298,-5.0725,2.8151,3.2492] ;
x_0 = x_0' ; 
disp("Initial gauss vale x_0 ") ;
disp(x_0);
s = eye(N) ; 

k = 1 ;
while true 
    for i=1:N 
        % step 2 ############################################################
        a = A ;
        b = B ;
        delta = 0.05 ;
        alfa_0 = a + (b-a)*rand  ; 
        
        % find mim point in along s1 direction ...........
        % Unidirectional Search: Bounding Phase..............
        
        f_0 = Objective_Fun(x_0 + alfa_0*s(:,i),N ) ;
        f_1 = Objective_Fun(x_0 + (alfa_0 + delta)*s(:,i)  ,N);
        f_2 = Objective_Fun(x_0 + (alfa_0 - delta)*s(:,i)  ,N);  
        
        
        if(f_1<f_0 && f_0<f_2)     % ......................check condition for detla
            delta = +delta ;
        end 
        if(f_2<f_0 && f_0<f_1)
            delta = -delta ;
        end 
        
        condition = true ;
        itb = 0 ; 
        while condition
            % step 3 
            alfa_1 = alfa_0 + 2.^itb * delta ;
            alfa_2 = alfa_1 + 2.^(itb+1) * delta ; 
            f1 = Objective_Fun(x_0 + alfa_1*s(:,i),N )  ;
            f2 = Objective_Fun(x_0 + alfa_2*s(:,i),N )  ;
          
            if(f2<f1)
                condition = true ;
                itb = itb + 1 ;
                
                alfa_0 = alfa_1 ;
            end     
            if(f2>f1)            % Termination condition
                condition = false ;
                itb = itb + 1 ;
            end     
        end 
        
        if(delta>0)
            a =alfa_0;
            b =alfa_2; 
        end 
        if(delta<0)
            a =alfa_2;
            b =alfa_0; 
        end 
        % Unidirectional Search: Interval Halving..............
        % step 1 
        e = 0.001 ;      % accuracy
        L = b - a ;      % Interval length
        alfa_m = (a+b)/2 ;   
        f_alfa_m = Objective_Fun(x_0 + alfa_m*s(:,i),N ) ;     % calling objective function
        
        %step 2 
        alfa1 = a+L/4 ;
        alfa2 = b-L/4 ;
        f_alfa1 = Objective_Fun(x_0 + alfa1*s(:,i),N ) ;
        f_alfa2 = Objective_Fun(x_0 + alfa2*s(:,i),N ) ;
        
        itE = 0 ;        % iteration in interval halving
        while(abs(L)>e)  
            % step 3
            if f_alfa1<f_alfa_m 
                b = alfa_m ;
                alfa_m = alfa1 ;
                
            elseif f_alfa2<f_alfa_m           
                 a = alfa_m ;
                alfa_m = alfa2 ;
            else 
                a=alfa1;
                b=alfa2;
            end
            itE = itE + 1 ;
            %fa = Objective_Fun(x_0 + a*s1,N);
            %fb = Objective_Fun(b);
            %fprintf(out,"%d %.5f\t  %.5f\t %.5f\t %.5f\n",itE,a,b,fa,fb);
            % step 4 
            % updating values 
            L = b-a ;
            alfa1 = a+L/4 ;
            alfa2 = b-L/4 ;
            f_alfa_m = Objective_Fun(x_0 + alfa_m*s(:,i),N );     % objective function calling
            f_alfa1 = Objective_Fun(x_0 + alfa1*s(:,i),N );
            f_alfa2 = Objective_Fun(x_0 + alfa2*s(:,i),N ) ;    
        end  
        alfa = alfa_m  ;
        x_1 = x_0 + alfa*s(:,i) ;  %......................new point x_1  
        x_old(:,i) = x_1 ;
        x_0 = x_1 ;
    end     
    
    % one more unidirection########################################### 
    
    
    % step 2 ############################################################
        i=1 ;
        a = A;
        b = B ;
        delta = 0.05 ;
        alfa_0 = a + (b-a)*rand  ; 
        
        % find mim point in along s1 direction ...........
        % Unidirectional Search: Bounding Phase..............
        
        f_0 = Objective_Fun(x_0 + alfa_0*s(:,i),N ) ;
        f_1 = Objective_Fun(x_0 + (alfa_0 + delta)*s(:,i)  ,N);
        f_2 = Objective_Fun(x_0 + (alfa_0 - delta)*s(:,i)  ,N);  
        
        
        if(f_1<f_0 && f_0<f_2)     % ......................check condition for detla
            delta = +delta ;
        end 
        if(f_2<f_0 && f_0<f_1)
            delta = -delta ;
        end 
        
        condition = true ;
        itb = 0 ; 
        while condition
            % step 3 
            alfa_1 = alfa_0 + 2.^itb * delta ;
            alfa_2 = alfa_1 + 2.^(itb+1) * delta ; 
            f1 = Objective_Fun(x_0 + alfa_1*s(:,i),N )  ;
            f2 = Objective_Fun(x_0 + alfa_2*s(:,i),N )  ;
          
            if(f2<f1)
                condition = true ;
                itb = itb + 1 ;
                
                alfa_0 = alfa_1 ;
            end     
            if(f2>f1)            % Termination condition
                condition = false ;
                itb = itb + 1 ;
            end     
        end 
        
        if(delta>0)
            a =alfa_0;
            b =alfa_2; 
        end 
        if(delta<0)
            a =alfa_2;
            b =alfa_0; 
        end 
        % Unidirectional Search: Interval Halving..............
        % step 1 
        e = 0.001 ;      % accuracy
        L = b - a ;      % Interval length
        alfa_m = (a+b)/2 ;   
        f_alfa_m = Objective_Fun(x_0 + alfa_m*s(:,i),N ) ;     % calling objective function
        
        %step 2 
        alfa1 = a+L/4 ;
        alfa2 = b-L/4 ;
        f_alfa1 = Objective_Fun(x_0 + alfa1*s(:,i),N ) ;
        f_alfa2 = Objective_Fun(x_0 + alfa2*s(:,i),N ) ;
        
        itE = 0 ;        % iteration in interval halving
        while(abs(L)>e)  
            % step 3
            if f_alfa1<f_alfa_m 
                b = alfa_m ;
                alfa_m = alfa1 ;
                
            elseif f_alfa2<f_alfa_m           
                 a = alfa_m ;
                alfa_m = alfa2 ;
            else 
                a=alfa1;
                b=alfa2;
            end
            itE = itE + 1 ;
            %fa = Objective_Fun(x_0 + a*s1,N);
            %fb = Objective_Fun(b);
            %fprintf(out,"%d %.5f\t  %.5f\t %.5f\t %.5f\n",itE,a,b,fa,fb);
            % step 4 
            % updating values 
            L = b-a ;
            alfa1 = a+L/4 ;
            alfa2 = b-L/4 ;
            f_alfa_m = Objective_Fun(x_0 + alfa_m*s(:,i),N );     % objective function calling
            f_alfa1 = Objective_Fun(x_0 + alfa1*s(:,i),N );
            f_alfa2 = Objective_Fun(x_0 + alfa2*s(:,i),N ) ;    
        end  
        alfa = alfa_m  ;
        x_1 = x_0 + alfa*s(:,i) ;  %......................new point x_1  
        x_old(:,N+1) = x_1 ;
        x_0 = x_1 ;
     
% ################################# find new direction  

    d = x_old(:,1) - x_old(:,N+1) ; 
    dn = norm(d) ; 
    if(dn<eps)
          break ;
    end    
    s_old = s ;
    for i=2:N
        s(:,i) = s_old(:,i-1) ;   
    end 
    s(:,1) = d/dn ; 
    k = 1+ k ;
end
fop = Objective_Fun(x_1,N);
disp("optimal point ");
disp(x_1);
disp("optimal value");
disp(fop);