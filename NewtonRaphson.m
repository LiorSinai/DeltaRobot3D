%%% 26 November 2018
%%% Lior Sinai
%%% Dynamics Problem Set 2
%%% Newton-Raphson method

function [qi,iter]=NewtonRaphson(symbols,Phi,Phi_q,qi,ti,epsilon)
    eqError=1;
    iter=0;
    maxIter=20;
    while(abs(eqError)>epsilon)
        iter=iter+1;
        Phi_i=subs(Phi,symbols,{[qi' ti]});
        eqError=double(norm(Phi_i));
        if (abs(eqError)<epsilon) 
            break; 
        end
        Phi_qi=subs(Phi_q,symbols,{[qi' ti]});
        qi=double((qi-Phi_qi\Phi_i)); 
        
        if iter>maxIter
            text(-2,-2,sprintf('Max iterations of %d exceeded',maxIter),'Color',[1 0 0])
            error('Max iterations of %d exceeded',maxIter);
        end
    end

end