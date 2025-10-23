function [c,ceq] = nlcon(betagammahat)
    c=[];
    ceq = norm(betagammahat(:))-1; 
end