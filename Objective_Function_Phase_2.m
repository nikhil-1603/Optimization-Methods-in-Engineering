% Q1 sumsquares function..........
function fun_val = Objective_Fun(x,N) 
    fun_val = 0 ;
    for i=1:N
    fun_val = fun_val +  i*(x(i))^2 ;
    end    
end  

% Q2 Rosenbrock function
function fun_val = Objective_Fun(x,N) 
    fun_val = 0 ;
    for i=1:N-1
    fun_val = fun_val +  100 * (x(i+1) - x(i)^2)^2 + (x(i) - 1)^2; 
    end    
end  

%Q3 Dixon Price function
function fun_val = Objective_Fun(x,N) 
    sum = 0 ;
    for i=2:N
    sum = sum +  i * (2*x(i)^2 - x(i-1))^2;
    end  
    fun_val = ( x(1) - 1)^2 + sum ;
end 



%Q4 Trid function 
function fun_val = Objective_Fun(x,N) 
    sum1 = 0 ;
    for i=2:N
        sum1 = sum1 + ( x(i) -1 )^2 ;
    end  
    sum2 = 0 ;
    for i=2:N
        sum2 = sum2 + x(i) * x(i-1) ;
    end    
    fun_val = sum1 - sum2 ;
end  



%Q5 Zakharov  function
function fun_val = Objective_Fun(x,N) 
    sum1 = 0 ;
    for i=1:N
        sum1 = sum1 + x(i)^2 ;
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
