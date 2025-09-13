clear all ;
clc;

function min_valu = min_fun(x)
    min_valu = -(  (sin(2*pi*x(1)))^3 * sin(2*pi*x(2))  ) / ( x(1)^3 *( x(1)+x(2) ) ) ; 
end     

function fun_val = Objective_Fun(x,R)     %Penalty function
    fun = -(  (sin(2*pi*x(1)))^3 * sin(2*pi*x(2))  ) / ( x(1)^3 *( x(1)+x(2) ) ) ;
    g1 = -x(1)^2 + x(2) -1 ;
    g2 = -1 + x(1) -(x(2) -4 )^2 ; 
    if(g1>0)
        g1 = 0 ;
    end 
    if(g2>0)
        g2 = 0 ;
    end  
    fun_val = fun + R*(g1*g1 + g2*g2) ;
end    

out = fopen('opt_phase3_Q2_out.txt', 'w'); % Output file
% step 1 ##########################################################
N = 2;%input('Enter no of variable N = ') ;
A = 0;%input('Enter lower  limit of x is A = ') ;
B = 10;%input('Enter upper  limit of x is B = ') ;
R(1) = 0.1 ;
eps = 0.001 ;
for i=1:N 
  x_0(i) = A + (B-A)*rand  ;    %  Random inital guess
end 
x_0 = [0.9428 4.0546] ;
x_0 = x_0' ;       % Transpose
disp("Initial gauss vale x_0 ") ;
disp(x_0);
s = eye(N) ;       % Identity matrix (initial direction)
k = 0 ;
count = 0 ;
fprintf(out,"%d\t %d\t\t %f\t %f\t %f\n",k,count,x_0(1),x_0(2),-min_fun(x_0)) ;

Ep = 5 ;

while Ep > 0.001

while true 
    for i=1:N 
        s(:,i) = s(:,i)/norm(s(:,i)) ;      % Unit direction conversion
        for j=1:N
            if(s(j,i)>0)                    % alpha bounds
                a_rang(j) = (A - x_0(j,1))/s(j,i) ;
                b_rang(j) = (B - x_0(j,1))/s(j,i) ;
            elseif(s(j,i)<0) 
                 b_rang(j) = (A - x_0(j,1))/s(j,i) ;
                 a_rang(j) = (B - x_0(j,1))/s(j,i) ;
            end    
        end 
        a = max(a_rang) ;
        b = min(b_rang) ;
        delta = 0.05 ;
        %alfa_0 = a + (b-a)*rand  ;
        alfa_0 =  0;
        
        % find min point along s direction ...........
        % Unidirectional Search Bounding Phase..............
        
        f_0 = Objective_Fun(x_0 + alfa_0*s(:,i),R(k+1)) ;
        f_1 = Objective_Fun(x_0 + (alfa_0 + delta)*s(:,i),R(k+1));
        f_2 = Objective_Fun(x_0 + (alfa_0 - delta)*s(:,i),R(k+1));  
        count = count + 3 ;
        
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
            f1 = Objective_Fun(x_0 + alfa_1*s(:,i),R(k+1) )  ;
            f2 = Objective_Fun(x_0 + alfa_2*s(:,i),R(k+1))  ;
          
            if(f2<f1)
                condition = true ;
                itb = itb + 1 ;
                
                alfa_0 = alfa_1 ;
            end     
            if(f2>f1)                    % Termination condition for Bounding Phase
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
        f_alfa_m = Objective_Fun(x_0 + alfa_m*s(:,i),R(k+1)) ;     
        
        %step 2 
        alfa1 = a+L/4 ;
        alfa2 = b-L/4 ;
        f_alfa1 = Objective_Fun(x_0 + alfa1*s(:,i),R(k+1)) ;
        f_alfa2 = Objective_Fun(x_0 + alfa2*s(:,i),R(k+1)) ;
        count = count + 3 ;
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
            % step 4 
            % updating values 
            L = b-a ;
            alfa1 = a+L/4 ;
            alfa2 = b-L/4 ;
            f_alfa_m = Objective_Fun(x_0 + alfa_m*s(:,i),R(k+1));   
            f_alfa1 = Objective_Fun(x_0 + alfa1*s(:,i),R(k+1));
            f_alfa2 = Objective_Fun(x_0 + alfa2*s(:,i),R(k+1)) ;  
            count = count + 3 ;
        end  
        alfa = alfa_m  ;
        x_1 = x_0 + alfa*s(:,i) ;  %......................new point  
        x_old(:,i) = x_1 ;
        x_0 = x_1 ;
    end     
    
    %   One more unidirection Search............... 
    
    
    % step 2 ############################################################
        i=1 ; 
        s(:,i) = s(:,i)/norm(s(:,i)) ;
        for j=1:N
            if(s(j,i)>0)
                a_rang(j) = (A - x_0(j,1))/s(j,i) ;
                b_rang(j) = (B - x_0(j,1))/s(j,i) ;
            elseif(s(j,i)<0) 
                 b_rang(j) = (A - x_0(j,1))/s(j,i) ;
                 a_rang(j) = (B - x_0(j,1))/s(j,i) ;
            end    
        end 
        a = max(a_rang) ;
        b = min(b_rang) ;
        delta = 0.05 ;
        %alfa_0 = a + (b-a)*rand  ;
        alfa_0 =  0 ;
        
        % find mim point in along s1 direction ...........
        % Unidirectional Search: Bounding Phase..............
        
        f_0 = Objective_Fun(x_0 + alfa_0*s(:,i),R(k+1)) ;
        f_1 = Objective_Fun(x_0 + (alfa_0 + delta)*s(:,i),R(k+1));
        f_2 = Objective_Fun(x_0 + (alfa_0 - delta)*s(:,i),R(k+1));  
        count = count + 3 ;
        
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
            f1 = Objective_Fun(x_0 + alfa_1*s(:,i),R(k+1))  ;
            f2 = Objective_Fun(x_0 + alfa_2*s(:,i),R(k+1))  ;
          
            if(f2<f1)
                condition = true ;
                itb = itb + 1 ;
                
                alfa_0 = alfa_1 ;
            end     
            if(f2>f1)              % Termination condition
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
        f_alfa_m = Objective_Fun(x_0 + alfa_m*s(:,i),R(k+1)) ;    
        
        %step 2 
        alfa1 = a+L/4 ;
        alfa2 = b-L/4 ;
        f_alfa1 = Objective_Fun(x_0 + alfa1*s(:,i),R(k+1)) ;
        f_alfa2 = Objective_Fun(x_0 + alfa2*s(:,i),R(k+1)) ;
        count = count + 3 ;
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
            % step 4 
            % updating values 
            L = b-a ;
            alfa1 = a+L/4 ;
            alfa2 = b-L/4 ;
            f_alfa_m = Objective_Fun(x_0 + alfa_m*s(:,i),R(k+1));    
            f_alfa1 = Objective_Fun(x_0 + alfa1*s(:,i),R(k+1));
            f_alfa2 = Objective_Fun(x_0 + alfa2*s(:,i),R(k+1)) ;    
        end  
        alfa = alfa_m  ;
        x_1 = x_0 + alfa*s(:,i) ;  %......................optimal point  
        x_old(:,N+1) = x_1 ;
        x_0 = x_1 ;
    
     
% ################################# 

    d = x_old(:,N+1) - x_old(:,1) ; 
    dn = norm(d) ;  
    
    
    

    if(dn<eps)           % Termination for Conjugate direction method
          break ;
    end    
    s_old = s ;
    for i=2:N
        s(:,i) = s_old(:,i-1) ;   %  Direction update
    end 
    s(:,1) = d/dn ; 
end 
k = k+1 ;
fprintf(out,"%d\t %d\t %f\t %f\t %f\n",k,count,x_0(1),x_0(2),-min_fun(x_0)) ;
R(k+1)= 10*R(k) ;
Xc(:,k) = x_0  ;
disp(Xc(:,k));
disp(-min_fun(Xc(:,k))) ;
%disp(Objective_Fun(Xc(:,k),R(k))) ;
if(k~=1)
    Ep = abs( Objective_Fun(Xc(:,k),R(k)) - Objective_Fun(Xc(:,k-1),R(k-1)) );
end 
if k>=20 
    break 
end 

end 

