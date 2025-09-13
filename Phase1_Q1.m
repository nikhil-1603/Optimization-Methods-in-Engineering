clear all;
clc;

%......... input( step 1 )
a  = input('Enter a = lower limit of x = ');
b  = input('Enter b = upper limit of x = ');
x_0 = input('Enter x_0 = initial guess value = ');
delta = input('Enter value of delta = '); 
%;  x_0 = a + (b-a)*rand      % Initial value 
out = fopen('Q1out.txt', 'w'); % Output file
fprintf(out,"Bounding phase method \n");
fprintf(out,"it  x_1\t\t  x_2\t   f(x_1)\t\t   f(x_2)\n");


function fun_val = Objective_Fun(x)
    fun_val = (2*x - 5).^4 - (x.^2 - 1).^3 ;  %  objective function f(x)
    
end

%step 2

f_0 = Objective_Fun(x_0) ;
f_1 = Objective_Fun(x_0+delta);
f_2 = Objective_Fun(x_0-delta); 

% ......................check condition for detla
if(f_1>f_0 && f_0>f_2)
    delta = +delta ;
end 
if(f_2>f_0 && f_0>f_1)
    delta = -delta ;
end 
 
condition = true ;
itb = 0 ;                         % iteration in bounding phase
fun_evaluation_B = 3 ;
while condition
    % step 3 
    x_1 = x_0 + 2.^itb * delta ;
    x_2 = x_1 + 2.^(itb+1) * delta ; 
    f1 = Objective_Fun(x_1) ;
    f2 = Objective_Fun(x_2) ;
  
    if(f2>f1)
        condition = true ;
        itb = itb + 1 ;
        fprintf(out,"%d %.3f\t %.3f\t %.3f\t %.3f \n",itb,x_1,x_2,f1,f2);
        x_0 = x_1 ;
    end     
    if(f2<f1)            % Termination condition
        condition = false ;
        itb = itb + 1 ;
        fprintf(out,"%d %.3f\t %.3f\t %.3f\t %.3f \n",itb,x_1,x_2,f1,f2);
    end 
    
    % step 4
    
    fun_evaluation_B = fun_evaluation_B + 1 ;    % no. of function evaluations
end 

if(delta>0)
    a =x_0;
    b =x_2; 
end 
if(delta<0)
    a =x_2;
    b =x_0; 
end 

% prtining the values in file
fprintf(out,"*************************************\n");
fprintf(out,"After bounding phase method \n");
fprintf(out,"Number of iteration in bounding phase method =  %d\n",itb);  % printing no. of iterations
fprintf(out,"a = %.3f\n",a);    
fprintf(out,"b = %.3f\n",b);
fprintf(out,"The intraval is ( %.3f , %.3f )\n",a,b);
fprintf(out,"Number of function evaluation in Bounding Phase Method = %d\n",fun_evaluation_B);


%.................x.............complete bounding phase method........x.....
%.................Interval Halving Method start..................
fprintf(out,"***********************************************\n");
fprintf(out,"Interval Halving Method start\n");
fprintf(out,"it\t  a\t\t\t b\t\t\t f(a)\t f(b)\n");


% step 1 
e = 0.001 ;      % accuracy
L = b - a ;      % Interval length
xm = (a+b)/2 ;   
fxm = Objective_Fun(xm) ;     % calling objective function

%step 2 
x1 = a+L/4 ;
x2 = b-L/4 ;
fx1 = Objective_Fun(x1) ;
fx2 = Objective_Fun(x2) ;

fun_evaluation_E = 0 ; 
itE = 0 ;        % iteration in interval halving
while(abs(L)>e)  
    % step 3
    if fx1>fxm 
        b = xm ;
        xm = x1 ;
        
    elseif fx2>fxm           
         a = xm ;
        xm = x2 ;
    else 
        a=x1;
        b=x2;
    end
    itE = itE + 1 ;
    fa = Objective_Fun(a);
    fb = Objective_Fun(b);
    fprintf(out,"%d %.5f\t  %.5f\t %.5f\t %.5f\n",itE,a,b,fa,fb);
    % step 4 
    % updating values 
    L = b-a ;
    x1 = a+L/4 ;
    x2 = b-L/4 ;
    fxm = Objective_Fun(xm);     % objective function calling
    fx1 = Objective_Fun(x1);
    fx2 = Objective_Fun(x2) ; 
    
    fun_evaluation_E = fun_evaluation_E + 2 ;   % no. of function evaluations
end
%  priting the values in file
fprintf(out,"***********************************************\n");
fprintf(out,"After Interval Halving Method \n");
fprintf(out,"Number of iteration in Interval Halving method =  %d\n",itE);
fprintf(out,'value of x = %.5f\n',x2);
fprintf(out,"maximum value of function at 'x' is %.5f\n",fx2);
fprintf(out,"Number of function evaluation in Interval Halving Method = %d\n",fun_evaluation_E);
