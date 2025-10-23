function result = getBaseFunVal(vec,i,k,U)
    result = [];
    for j = 1:length(vec)
        u = vec(j);
        Val = 0;
        val1 = 0; val2 = 0;
         
        if (k == 0)
            if (u>=U(i +1))&&(u<U(i+1 +1))  
                Val = 1;
                Bf = Val;
                
            else
                Bf = Val;
            end
        end
         
        if (k>0)
            if  (u<U(i +1))||(u>=U(i+k+1 +1))
                Bf = Val;
            else
                alpha = 0;
                beta = 0;
                dTemp = 0;
                
                dTemp = U(i+k +1) - U(i+ 1);
                if dTemp == 0
                    alpha = 0;        
                else
                    alpha = (u - U(i+ 1))/dTemp;
                end
                
                dTemp = U(i+k+1 +1) - U(i+1 +1);
                if dTemp == 0
                    beta = 0;           
                else
                    beta = (U(i+k+1 +1) - u)/dTemp;
                end
                
                val1 = alpha*getBaseFunVal(u,i,k-1,U);  
                val2 = beta*getBaseFunVal(u,i+1,k-1,U); 
                Val = val1 + val2;
                Bf = Val;
            end
        end
        result = [result,Bf];
    end
