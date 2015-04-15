
function [u, l2_err] = SB_SPA1D(R,f, mu, lambda, gamma, nInner, nBreg,m_spa,UT)
  
   [rows,cols] = size(f(:));
    
        
    % Reserve memory for the auxillary variables
    N = rows;
    C_spa = ones(m_spa+1,1);
    for j = 1:m_spa+1
        for js =1:m_spa+1
            if js ~= j
                C_spa(j) = C_spa(j)/(j-js);
            end
        end
    end
    L_SPA = zeros(N);
    m2 = floor((m_spa+1)/2);
    q_norm = sum(C_spa(1:m2,1));
    C_spa = C_spa./q_norm ;
    for j=1:N
        inds = j-m2:j+m_spa-m2;
        inds(find(inds<1))=inds(find(inds<1))+N;
        inds(find(inds>N))=inds(find(inds>N))-N;
        L_SPA(j,inds) = C_spa;
    end

   % L_SPAT = transpose(L_SPA);
    
    f0 = f;
    u = zeros(rows,cols);
    x = zeros(rows,cols);
    bx = zeros(rows,cols);
    

     % Build Kernels
    scale = sqrt(rows*cols);
    murf = ifft(mu*(conj(R).*f))*scale;
    
    uker = zeros(rows,cols);
    uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
    uker = mu*(conj(R).*R)+lambda*fft(uker).^m_spa+gamma;

    l2_err = [];

    for outer = 1:nBreg;
        for inner = 1:nInner;
             % update u  
            rhs = murf+lambda*DR(x-bx,L_SPA)+gamma*u;
            u = ifft(fft(rhs)./uker);

            % update x and y
            dx = DR(u,L_SPAT);
            [x] = shrink( dx+bx,1/lambda);

            % update bregman parameters
            bx = bx+dx-x;
            
        end
        l2_err = [l2_err;sqrt(sum(sum(real(UT-u).^2)))/rows^2];
        f = f+f0-R.*fft(u)/scale;
        murf = ifft(mu*R.*f)*scale;
    end
    figure, plot(l2_err); drawnow
return;


function d = DR(u,M)
d = u*M;
return



function [xs] = shrink(x,lambda)

s = sqrt(x.*conj(x));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;

return;