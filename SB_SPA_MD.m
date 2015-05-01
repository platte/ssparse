
function u = SB_SPA_MD(R,f, mu, lambda, gamma, nInner, nBreg,m_spa)
    [rows,cols] = size(f);
    
        
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
    L_spa = zeros(N);
    m2 = floor((m_spa+1)/2);
    q_norm = sum(C_spa(1:m2,1));
    C_spa = C_spa./q_norm ;
    for j=1:N
        inds = j-m2:j+m_spa-m2;
        inds(find(inds<1))=inds(find(inds<1))+N;
        inds(find(inds>N))=inds(find(inds>N))-N;
        L_SPA(j,inds) = C_spa;
    end

    L_SPAT = transpose(L_SPA);
    
    f0 = f;
    u = zeros(rows,cols);
    x = zeros(rows,cols);
    y = zeros(rows,cols);
    bx = zeros(rows,cols);
    by = zeros(rows,cols);
    
    l2_err = [];

     % Build Kernels
    scale = sqrt(rows*cols);
    murf = ifft2(mu*(conj(R).*f))*scale;
    
    uker = zeros(rows,cols);
    uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
    uker = mu*(conj(R).*R)+lambda*fft2(uker).^m_spa+gamma;

    h_fig = figure;
    for outer = 1:nBreg;
        for inner = 1:nInner;
             % update u  
            up = u;
            rhs = murf+lambda*DR(x-bx,L_SPA)+lambda*DL(y-by,L_SPAT)+gamma*u;
            u = ifft2(fft2(rhs)./uker);

            % update x and y
            dx = DR(u,L_SPAT);
            dy  =DL(u,L_SPA);
            [x,y] = shrink2( dx+bx, dy+by,1/lambda);

            % update bregman parameters
            bx = bx+dx-x;
            by = by+dy-y;
        end
        l2_err = [l2_err;sqrt(sum(sum(real(up-u).^2)))/rows^2];
        f = f+f0-R.*fft2(u)/scale;
        murf = ifft2(mu*R.*f)*scale;
        figure(h_fig), semilogy(l2_err); drawnow
    end
    close(h_fig);
return;


function d = DL(u,M)
d = M*u;
return

function d = DR(u,M)
d = u*M;
return



function [xs,ys] = shrink2(x,y,lambda)

s = sqrt(x.*conj(x)+y.*conj(y));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;

return;


