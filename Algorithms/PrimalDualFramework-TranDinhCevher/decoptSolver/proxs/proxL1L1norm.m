function x = proxL1L1norm(z, w, beta, gamma)

v = -(1/gamma)*z;

rho = 1/gamma;

rhow = rho*w;                                                                   

n = length(w);

x = zeros(n,1);
                                                                                
w_pos = (w >= 0);


% ************************************************************  
% Components for which w_i >= 0

case1 = logical((w_pos) .* (v < -rhow -beta -1)); 
x(case1) = ( -beta - 1 - v(case1) )/rho; 

case2 = logical(w_pos .* (-rhow - beta - 1 <= v) .* (v <= -rhow + beta - 1));
x(case2) =  w(case2);

case3 = logical(w_pos .* (-rhow + beta - 1 < v) .* (v < beta - 1));
x(case3) = (beta - 1 - v(case3))/rho;

case4 = logical(w_pos .* (beta - 1 <= v) .* (v <= beta + 1));
x(case4) = 0;
    
case5 = logical(w_pos .* (v > beta + 1));
x(case5) = (beta + 1 - v(case5))/rho;
% ************************************************************  


% ************************************************************  
% Components for which w_i < 0

case1 = logical(~w_pos .* (v < -beta -1));
x(case1) = (-beta -1 - v(case1))/rho;

case2 = logical(~w_pos .* (-beta - 1 <= v) .* (v <= -beta + 1));
x(case2) = 0;

case3 = logical(~w_pos .* (-beta + 1 < v) .* (v < -rhow - beta + 1));
x(case3) = (-beta + 1 - v(case3))/rho;

case4 = logical(~w_pos .* (-rhow - beta + 1 <= v) .* (v <= -rhow + beta + 1));
x(case4) = w(case4);

case5 = logical(~w_pos .* (v > -rhow + beta + 1));
x(case5) = (beta + 1 - v(case5))/rho;
% ************************************************************  


% for i = 1 : n                                                                
%                                                                                 
%     if w(i) >= 0                                                                
%         if v(i) < -rhow(i) - beta - 1                                           
%             x(i) = (-beta -1 -v(i))/rho;                                        
%                                                                                 
%         elseif (-rhow(i) - beta - 1 <= v(i)) && (v(i) <= -rhow(i) + beta - 1)   
%             x(i) = w(i);                                                        
%                                                                                 
%         elseif (-rhow(i) + beta - 1 < v(i)) && (v(i) < beta - 1)                
%             x(i) = (beta - 1 -v(i))/rho;                                        
%                                                                                 
%         elseif (beta - 1 <= v(i)) && (v(i) <= beta + 1)                         
%             x(i) = 0;                                                           
%                                                                                 
%         else                                                                    
%             x(i) = (beta + 1 - v(i))/rho;                                       
%         end                                                                     
%                                                                                 
%     else                                                                        
%                                                                                 
%         if v(i) < -beta -1                                                      
%             x(i) = (-beta -1 - v(i))/rho;                                       
%                                                                                 
%         elseif (-beta - 1 <= v(i) )&& (v(i) <= -beta + 1)                       
%             x(i) = 0;                                                           
%                                                                                 
%         elseif (-beta + 1 < v(i)) && (v(i) < -rhow(i) - beta + 1)               
%             x(i) = (-beta + 1 - v(i))/rho;                                      
%                                                                                 
%         elseif (-rhow(i) - beta + 1 <= v(i)) && (v(i) <= -rhow(i) + beta + 1)   
%             x(i) = w(i);                                                        
%                                                                                 
%         else                                                                    
%             x(i) = (beta + 1 - v(i))/rho;                                       
%         end                                                                     
%                                                                                 
%     end                                                                         
% end                                                                             
