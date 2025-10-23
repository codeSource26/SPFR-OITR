% 产生数据
function[y,beta,z,gammat,x,A] = simulation(n,tp,sigma)
    % scalar covariates
    z1 = -1 + (1-(-1)).*rand(n,1); 
    z2 = -1 + (1-(-1)).*rand(n,1);
    z = [z1,z2]; 
    % functional covariates
    [x] = brownien(n,tp);
    % firstpart
    alpha = [0.3,0.4];
    t = linspace(0,1,tp);
    thetat = (0.5*t-1).^3 + 1;
    firstpart = (z*alpha' + x*thetat'*0.01);
    % secondpart
    beta = [1/sqrt(4),-1/sqrt(4)];
    gammat = sin(pi*t/2)/sqrt(2) + sin(3*pi*t/2)/sqrt(2);
    secondpart = 2*(z*beta' + x*gammat'*0.01).^3 + 0.6;
    % A
    A = binornd(1,0.5,n,1);
    A(A==0) = -1;
    epsnorm = normrnd(0,sigma,[n,1]);
    % total
    y = firstpart + secondpart.*A + epsnorm;
end









