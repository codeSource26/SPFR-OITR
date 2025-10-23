function[A] = Amatrix(N)
    A = zeros(N-1,N);
    for i = 1:(N-1)
        for j = 1:(N)
            if i == j
                A(i,j) = 1;
            elseif j == i+1
                A(i,j) = -1;
            else
                A(i,j) = 0;
            end
        end
    end
end
