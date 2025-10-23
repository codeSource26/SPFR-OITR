function N = Nmatrix(uhat,nknots,norder,nbasis)
    knots = quantile(uhat,linspace(0,1,nknots));
    knots = [knots(1)*ones(1,norder+1),knots(2:(length(knots)-1)),knots(length(knots))*ones(1,norder+1)];
    
    N = zeros(length(uhat),nbasis);
    for j = 1:nbasis
        Nj = getBaseFunVal(uhat,j,norder-1,knots); 
        N(:,j) = Nj';
        
        if j==nbasis
            [maxValue, linearIndex] = max(uhat);
            N(linearIndex,j) = 1;
        end
    end
end