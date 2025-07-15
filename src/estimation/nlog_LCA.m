function [lnL,g,h] = nlog_LCA(Params,Data,tau)
   
    [m,N] = size(Data); 
    leng = m - 1;
    
    xB = Data(1:end-1, :);
    xF = Data(2:end, :);
    
    xi_squared = Params(3);
    
    I = Params(4:end);
    barI = sum(I)/N;
    delI = I - barI;
    
    kbdiff = Params(1);
    kbNsum = Params(2);
    
    f0form = @(p) (1-exp(-p*tau))/p;
    b_5 = xi_squared / 2 * f0form(2*kbdiff);
    b_7 = xi_squared / 2 * f0form(2*kbNsum);
    
    invOmega = ones(N) * (1/(N*b_7)-1/(N*b_5)) + eye(N) / b_5;
    logdetO = (N-1)*log(b_5) + log(b_7);
    
    presum = sum(xB,2);
    y = (xB .* exp(-kbdiff*tau) + barI * f0form(kbNsum) + delI * f0form(kbdiff) + (exp(-kbNsum*tau)-exp(-kbdiff*tau))/N .* presum ) - xF ;
    
    A = reshape(y',1,N,leng);
    B = permute(A,[2 1 3]);
    C = pagemtimes(A,pagemtimes (invOmega,B));
    
    lnL = 0.25 * sum(C) + logdetO * leng * 0.5 ;

    if nargout > 1 % gradient required
        f1form = @(p) (exp(-p * tau)*(1 + p * tau)-1)/(p*p);
        db_5xi_squared = 0.5 * f0form(2*kbdiff);
        db_7xi_squared = 0.5 * f0form(2*kbNsum);

        db_5kbdiff = (xi_squared) * f1form(2*kbdiff);
        db_5kbNsum = 0;

        db_7kbdiff = 0;
        db_7kbNsum = (xi_squared) * f1form(2*kbNsum);
        db_5 = [db_5kbdiff db_5kbNsum db_5xi_squared ];
        db_7 = [db_7kbdiff db_7kbNsum db_7xi_squared ];

        dykbdiff = (-tau * xB .* exp(-kbdiff*tau) + delI * f1form(kbdiff) + (tau*exp(-kbdiff*tau))/N .* presum );
        dykbNsum = (0 * xB + barI * f1form(kbNsum) - (tau*exp(-kbNsum*tau))/N .* presum );
        dyxi = 0 * xB;
        dy = cat(3,dykbdiff,dykbNsum,dyxi);
        dA_i = @(j) reshape(dy(:,:,j)',1,N,leng);
        dB_i = @(j) reshape(dy(:,:,j)',N,1,leng);

        dyIj = ones(leng, N) * (1/N * f0form(kbNsum) - 1/N * f0form(kbdiff) );
        dyIi = ones(leng, 1) * (1/N * f0form(kbNsum) + (1-1/N) * f0form(kbdiff));
        dI = NaN(1,N);
        dlogdetO = (N-1)/b_5.*db_5 + 1/b_7 .* db_7;
        dinvOmega = kron((-1/(N*b_7^2) * db_7 + 1/(N*b_5^2) * db_5), ones(N)) - kron(1/b_5^2 * db_5  , eye(N));
        dinvOmega_j = @(j) dinvOmega(:,N*(j-1)+1:N*j);

        for i = 1:N
            dyI = dyIj;
            dyI(:, i) = dyIi;
            dA_I = reshape(dyI',1,N,leng);

            dC_Iy = 2 * pagemtimes(dA_I,pagemtimes(invOmega,B)); 
            dI(i) = 0.25 * sum(dC_Iy);
        end
        dC = NaN(1,3);
        for j = 1:3
            dC(j) = sum(2*pagemtimes(dA_i(j),pagemtimes(invOmega,B)) + pagemtimes(A,pagemtimes(dinvOmega_j(j),B)));
        end
        dparams = 0.25 * dC + dlogdetO * leng * 0.5;
        g = [dparams,dI];
        if nargout > 2
            f2form = @(p) -1*(exp(-p * tau)*(2 + 2 *  p * tau + p^2 * tau^2)-2)/(p^3);
            ddb_5kbdiffkbdiff = 2 * xi_squared * f2form(2*kbdiff);
            ddb_5kbdiffkbNsum = 0;
            ddb_5kbNsumkbNsum = 0;
            ddb_5kbdiffxi_squared = f1form(2*kbdiff);
            ddb_5kbNsumxi_squared = 0;
            ddb_5xi_squaredxi_squared = 0; 

            ddb_5 = [ddb_5kbdiffkbdiff,0,0;ddb_5kbdiffkbNsum,ddb_5kbNsumkbNsum,0;ddb_5kbdiffxi_squared,ddb_5kbNsumxi_squared,ddb_5xi_squaredxi_squared];
            ddb_5 = tril(ddb_5,-1)'+ddb_5;

            ddb_7kbdiffkbdiff = 0;
            ddb_7kbdiffkbNsum = 0;
            ddb_7kbNsumkbNsum = 2 * xi_squared * f2form(2*kbNsum);
            ddb_7kbdiffxi_squared = 0;
            ddb_7kbNsumxi_squared = f1form(2*kbNsum); 
            ddb_7xi_squaredxi_squared = 0; 

            ddb_7 = [ddb_7kbdiffkbdiff,0,0;ddb_7kbdiffkbNsum,ddb_7kbNsumkbNsum,0;ddb_7kbdiffxi_squared,ddb_7kbNsumxi_squared,ddb_7xi_squaredxi_squared];
            ddb_7 = tril(ddb_7,-1)'+ddb_7;

            db_5ij = db_5'*db_5;
            db_7ij = db_7'*db_7;

            h = zeros(N+3);
            ddlogdetO = (N-1) * (-1/b_5^2 * db_5ij + 1/b_5 * ddb_5) + (-1/b_7^2 * db_7ij + 1/b_7 * ddb_7);

            ddinvOmega = kron(1/N*(2/(b_7^3) * db_7ij - 1/(b_7^2) * ddb_7) - 1/N*(2/(b_5^3) * db_5ij - 1/(b_5^2) * ddb_5),ones(N)    ) + kron( (2/(b_5^3) * db_5ij - 1/(b_5^2) * ddb_5),eye(N)   );
            ddinvOmega = reshape(ddinvOmega, [N, 3, N, 3]);
            ddinvOmega = permute(ddinvOmega, [1, 3, 2, 4]);
            ddinvOmega_ij = @(i,j) reshape(ddinvOmega(:,:,i,j), N, N);
            ddinvOmega2 = kron(1/N*(2/(b_7^3) * db_7ij - 1/(b_7^2) * ddb_7) - 1/N*(2/(b_5^3) * db_5ij - 1/(b_5^2) * ddb_5),ones(N)    ) + kron( (2/(b_5^3) * db_5ij - 1/(b_5^2) * ddb_5),eye(N)   );
            ddinvOmega2_ij = @(i,j) ddinvOmega2(N*(i-1)+1:N*i,N*(j-1)+1:N*j);
            
            ddykbdiff = (tau^2 * xB .* exp(-kbdiff*tau) + delI * f2form(kbdiff) - (tau^2*exp(-kbdiff*tau))/N .* presum );
            ddykbNsum = (0 * xB + barI * f2form(kbNsum) + (tau^2*exp(-kbNsum*tau))/N .* presum );
            ddyxi = 0 * xB;
            ddy = cat(3,ddykbdiff,ddykbNsum,ddyxi);
            ddyIjkbdiff = ones(leng, N) * ( -1/N * f1form(kbdiff) );
            ddyIjkbNsum = ones(leng, N) * (  1/N * f1form(kbNsum) );
            ddyIjxi = zeros(leng, N);
            ddyIj = cat(3,ddyIjkbdiff, ddyIjkbNsum, ddyIjxi);
            ddyIikbdiff = ones(leng, 1) * ( (1-1/N) * f1form(kbdiff) );
            ddyIikbNsum = ones(leng, 1) * (    1/N  * f1form(kbNsum) );
            ddyIixi = zeros(leng, 1);
            ddyIi = cat(3,ddyIikbdiff, ddyIikbNsum, ddyIixi);
            hAXA = NaN(3,3);
            hAXA2 = NaN(3,3);
            for r=1:3
                for s=1:3
                    if (r==s)
                        ddA = reshape(ddy(:,:,r)',1,N,leng);
                        ddB = reshape(ddy(:,:,r)',N,1,leng);
                    else
                        ddA = 0 * zeros(1,N,leng);
                        ddB = 0 * zeros(N,1,leng);
                    end                   
                    
                    dAr = dA_i(r);
                    dBr = dB_i(r);
                    hAXA(r,s) = sum(...
                        + pagemtimes(ddA,pagemtimes(invOmega,B)) ...
                        + pagemtimes(A,pagemtimes(invOmega,ddB)) ...
                        + pagemtimes(dA_i(r),pagemtimes(dinvOmega_j(s),B))...
                        + pagemtimes(A,pagemtimes(dinvOmega_j(s),dB_i(r)))...
                        + pagemtimes(dA_i(r),pagemtimes(invOmega,dB_i(s)))...
                        + pagemtimes(dA_i(s),pagemtimes(invOmega,dB_i(r)))...
                        + pagemtimes(dA_i(s),pagemtimes(dinvOmega_j(r),B))...
                        + pagemtimes(A,pagemtimes(dinvOmega_j(r),dB_i(s)))...
                        + pagemtimes(A,pagemtimes(ddinvOmega_ij(r,s),B))...
                        );
                 end
            end
            for p = 1:N
                dyIp = dyIj;
                dyIp(:,p) = dyIi;
                dA_Ip = reshape(dyIp',1,N,leng);
                dB_Ip = reshape(dyIp',N,1,leng);
                for q = p:N
                    dyIq = dyIj;
                    dyIq(:,q) = dyIi;
                    dA_Iq = reshape(dyIq',1,N,leng);
                    dB_Iq = reshape(dyIq',N,1,leng);
                    dC_Iy = pagemtimes(dA_Ip,pagemtimes(invOmega,dB_Iq)) + pagemtimes(dA_Iq,pagemtimes(invOmega,dB_Ip)); 
                    h(3+p,3+q) = 0.25 * sum(dC_Iy);
                end
            end
            for s = 1:3 
                for p = 1:N
                    dyIp = dyIj;
                    dyIp(:,p) = dyIi;
                    dA_Ip = reshape(dyIp',1,N,leng);
                    dB_Ip = reshape(dyIp',N,1,leng);
                    ddyIAp = ddyIj(:,:,s);
                    ddyIAp(:, p) = ddyIi(:,:,s);
                    ddyIAp = reshape(ddyIAp',1,N,leng);
                    offdiagonal = 0.25 * 2 * sum(...
                        + pagemtimes(ddyIAp,pagemtimes(invOmega,B)) ...
                        + pagemtimes(dA_i(s),pagemtimes(invOmega,dB_Ip)) ...
                        + pagemtimes(dA_Ip,pagemtimes(dinvOmega_j(s),B)) ...
                        );
                    h(s,p+3) = offdiagonal;
                end
            end
            h(1:3,1:3) = 0.25 * hAXA + ddlogdetO * leng * 0.5;
            h = triu(h)+triu(h,1)';
            
        end
    end
end