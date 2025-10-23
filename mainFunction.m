n = 1000;
sigma = 0.4;
M = 500; %repetition number
L = 500; %maximum number of iterations
p = 2; %dimension of beta
mMax = 5; %max value of m
tp = 100; %number of discrete points of x
norder = 4; %degree of B-spline
nknotsMax = 12; %maximum number of knots



% ======simulation start========================
zplot = linspace(-1,1,101);
Bspline_betaerror = zeros(M,p);
Bspline_betaAbserror = zeros(M,p);
Bspline_betaMse = zeros(M,p);
Bspline_psierror = zeros(1,M);
Bspline_zplotPre = zeros(M,length(zplot));
Bspline_pcd = zeros(1,M);
Bspline_valEst = zeros(1,M);
Bspline_valEstError = zeros(1,M);
Bspline_alist = zeros(1,M);
Bspline_regret = zeros(1,M);
Bspline_gammatMse = zeros(1,M);
Bspline_GammaCurveHat = zeros(M,tp);
betaHatStdHat = zeros(M,p);
IntervalCover = zeros(M,p);
mVec = zeros(1,M);
nknotsVec = zeros(1,M);


for k = 1:M
    % generate data
    [y1,beta,z,gammat,x,A] = simulation(n,tp,sigma);
    % yOriginal = y;
    w = zeros(n,1);
    w(A==1) = A(A==1)/(2*(sum(A(A==1))/n));
    w(A==-1) = A(A==-1)/(2*(-sum(A(A==-1))/n));
    y = y1.*w;

    BIC_m_nknots = zeros(mMax,nknotsMax);

    for nknots = 5:nknotsMax
        nbasis = nknots + norder - 2;
        for m = 1:mMax
            %pca
            [V,D] = eig(cov(x));
            v = zeros(tp,m);
            score = zeros(n,m);
            for j = 1:m
                v(:,j) = sqrt(tp)*V(:,end-j+1); 
                score(:,j) = x*v(:,j)*0.01;
            end
            Regressor = [z,score]; 

            %------------------------------Bspline method--------------------------------------
            %Step1 compute the initial value of beta
            betahat0 = lsqlin(Regressor,y);
            betahat0 = normalize(betahat0,'norm');
        
            Bspline_alphahat = zeros(L+1,nbasis);
            Bspline_betagammahat = zeros(L+1,length(beta)+m);
            Bspline_betagammahat(1,:) = betahat0;
        
        
            for l = 1:L
                %Step2 Quadratic programming, to solve alpha
                Ne = Nmatrix(Regressor*Bspline_betagammahat(l,:)',nknots,norder,nbasis);
                [alphahattemp,fval] = quadprog(Ne'*Ne,-Ne'*y,AmatrixBspline(nbasis),zeros(1,nbasis-1),[],[],[],[],[],optimoptions('quadprog','Display','off'));
                Bspline_alphahat(l+1,:) = alphahattemp;
        
                %Step3 Nonlinear least squares, to solve beta
                fun = @(betagammahat)Nmatrix(Regressor*betagammahat',nknots,norder,nbasis)*alphahattemp - y;
                betagammahattemp = lsqnonlin(fun,Bspline_betagammahat(l,:),[],[],[],[],[],[],@nlcon,optimoptions('lsqnonlin','Display','off'));
                Bspline_betagammahat(l+1,:) = betagammahattemp;
                Bspline_betagammahat(l+1,:)=Bspline_betagammahat(l+1,:)/sqrt(sum(Bspline_betagammahat(l+1,:).^2));
                
                %determine if convergence
                a = l; %record the number of times each simulation iteration converges
                if (norm(Bspline_betagammahat(l+1,:)-Bspline_betagammahat(l,:))<1e-3 && norm(Bspline_alphahat(l+1,:)-Bspline_alphahat(l,:))<1e-3)
                    break;
                end
            end
        
       
            psierror = sum((Nmatrix(Regressor*Bspline_betagammahat(a+1,:)',nknots,norder,nbasis)*Bspline_alphahat(a+1,:)' - (2*(z*beta' + x*gammat'*0.01).^3 + 0.6)).^2)/n; 
            BIC_m_nknots(m,nknots) = n*log(psierror) + (nknots + m + nknots)*log(n);

        end
    end

    [m,nknots] = find(BIC_m_nknots==min(min(BIC_m_nknots(:,5:end))));
    nbasis = nknots + norder - 2; 
    mVec(k) = m;
    nknotsVec(k) = nknots;


    Bspline_BetaGammaHat = zeros(M,p+m);
    Bspline_BetaGammaHat0 = zeros(M,p+m);
    Bspline_AlphaHat = zeros(M,nbasis);
    % pca
    [V,D] = eig(cov(x));
    v = zeros(tp,m);
    score = zeros(n,m);
    for j = 1:m
        v(:,j) = sqrt(tp)*V(:,end-j+1); 
        score(:,j) = x*v(:,j)*0.01;
    end

    Regressor = [z,score]; 

    %------------------------------estimate using optimal parameters--------------------------------------
    %Step1 compute the initial value of beta
    betahat0 = lsqlin(Regressor,y);
    betahat0 = normalize(betahat0,'norm');
    Bspline_BetaGammaHat0(k,:) = betahat0;

    Bspline_alphahat = zeros(L+1,nbasis);
    Bspline_betagammahat = zeros(L+1,length(beta)+m);
    Bspline_betagammahat(1,:) = betahat0;


    for l = 1:L
        %Step2 Quadratic programming, to solve alpha
        Ne = Nmatrix(Regressor*Bspline_betagammahat(l,:)',nknots,norder,nbasis);
        [alphahattemp,fval] = quadprog(Ne'*Ne,-Ne'*y,AmatrixBspline(nbasis),zeros(1,nbasis-1),[],[],[],[],[],optimoptions('quadprog','Display','off'));
        Bspline_alphahat(l+1,:) = alphahattemp;

        %Step3 Nonlinear least squares, to solve beta
        fun = @(betagammahat)Nmatrix(Regressor*betagammahat',nknots,norder,nbasis)*alphahattemp - y;
        betagammahattemp = lsqnonlin(fun,Bspline_betagammahat(l,:),[],[],[],[],[],[],@nlcon,optimoptions('lsqnonlin','Display','off'));
        Bspline_betagammahat(l+1,:) = betagammahattemp;
        Bspline_betagammahat(l+1,:)=Bspline_betagammahat(l+1,:)/sqrt(sum(Bspline_betagammahat(l+1,:).^2));
        
        %determine if convergence
        a = l; %record the number of times each simulation iteration converges
        if (norm(Bspline_betagammahat(l+1,:)-Bspline_betagammahat(l,:))<1e-3 && norm(Bspline_alphahat(l+1,:)-Bspline_alphahat(l,:))<1e-3)
            break;
        end
    end


    Bspline_alist(k) = a;
    Bspline_BetaGammaHat(k,:) = Bspline_betagammahat(a+1,:);
    Bspline_AlphaHat(k,:) = Bspline_alphahat(a+1,:);

    Bspline_betaerror(k,:) = Bspline_BetaGammaHat(k,1:p) - beta; 
    Bspline_betaAbserror(k,:) = abs(Bspline_BetaGammaHat(k,1:p) - beta); 
    Bspline_betaMse(k,:) = (Bspline_BetaGammaHat(k,1:p) - beta).^2; 
    Bspline_psierror(k) = sum((Nmatrix(Regressor*Bspline_betagammahat(a+1,:)',nknots,norder,nbasis)*Bspline_alphahat(a+1,:)' - (2*(z*beta' + x*gammat'*0.01).^3 + 0.6)).^2)/n; 
    
    % predict of psi on zplot
    uhat = Regressor*Bspline_betagammahat(a+1,:)';
    modified_uhat = uhat;
    if max(uhat) < 1.01
        [max_val, max_idx] = max(uhat);
        modified_uhat(max_idx) = 1.01;
    end
    if min(uhat) > -1.01
        [min_val, min_idx] = min(uhat);
        modified_uhat(min_idx) = -1.01;
    end
    knots = quantile(modified_uhat,linspace(0,1,nknots));
    knots = [knots(1)*ones(1,norder+1),knots(2:(length(knots)-1)),knots(length(knots))*ones(1,norder+1)];
    Ne = zeros(length(zplot),nbasis);
    for j = 1:nbasis
        Nej = getBaseFunVal(zplot,j,norder-1,knots); 
        Ne(:,j) = Nej';
    end
    Bspline_zplotPre(k,:) = Ne*Bspline_alphahat(a+1,:)';


    % MSE of gammat
    Bspline_gammatMse(k) = sum((Bspline_BetaGammaHat(k,(p+1):(p+m))*v' - gammat).^2) * 0.01; 
    % estimate of gammat
    Bspline_GammaCurveHat(k,:) = Bspline_BetaGammaHat(k,(p+1):(p+m))*v';

    % pcd of the estimated optimal itr
    itrEst = sign(Nmatrix(Regressor*Bspline_betagammahat(a+1,:)',nknots,norder,nbasis)*Bspline_alphahat(a+1,:)');
    itrTrue = sign(2*(z*beta' + x*gammat'*0.01).^3 + 0.6);
    Bspline_pcd(k) = 1 - sum(abs(itrTrue - itrEst))/(2*n);  

    % compute |E[Y({\hat d}^{opt})]-E[Y(d^{opt})]|
    Bspline_regret(k) = abs(sum((2*(z*beta' + x*gammat'*0.01).^3 + 0.6).*(sign(Nmatrix(Regressor*Bspline_betagammahat(a+1,:)',nknots,norder,nbasis)*Bspline_alphahat(a+1,:)') - sign(2*(z*beta' + x*gammat'*0.01).^3 + 0.6))) / n);



    k
end

%(1)bias of betahat
Bspline_betaerrorMean = mean(Bspline_betaerror,1);
Bspline_betaerrorMean
%(2)absolute bias of betahat
Bspline_betaAbserrorMean = mean(Bspline_betaAbserror,1);
Bspline_betaAbserrorMean
%(3)MSE of betahat
Bspline_betaMseMean = mean(Bspline_betaMse,1);
Bspline_betaMseMean
%(4)MSE of gammat
Bspline_gammatMseMean = mean(Bspline_gammatMse);
Bspline_gammatMseMean
%(5)MSE of psihat
Bspline_psierrorMean = mean(Bspline_psierror);
Bspline_psierrorMean
%(6)pcd
Bspline_pcdMean = mean(Bspline_pcd);
Bspline_pcdMean
%(7)regret
Bspline_regretMean = mean(Bspline_regret);
Bspline_regretMean


%summarize results================================================================================
[Bspline_betaerrorMean,Bspline_betaAbserrorMean,Bspline_betaMseMean,Bspline_gammatMseMean,Bspline_psierrorMean,Bspline_pcdMean,Bspline_regretMean]














