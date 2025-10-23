 function [x] = brownien(n,tp)
    x = zeros(n,tp);
    for i = 1:n
        x(i,:) = cumsum((1/sqrt(tp))*randn(1,tp)); 
    end
end


