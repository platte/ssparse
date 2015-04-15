
close all
clear all
clc

N = 128; % The image will be NxN
sparsity = .2; % use only 50% on the K-Space data for CS
Tstar = 4*pi;
tvar = 1:5*pi:10*pi;
R2D2 = @(x,y,t) -t/Tstar;
FreqF = @(x,y,t) .25*(exp(-12*((y-.6).^2+x.^2))+exp(-12*((y+.6).^2+x.^2)))*t;
Window_Patch = 3;
n_inner = 10;
n_outer = 200;
n_inner_spa = 10;
n_outer_spa = 10;
n_inner_spa3 = 10;
n_outer_spa3 = 10;

mu = 1;
mu_spa = 1;
lambda = 100;

gamma_tv = 1;



% build an image of a SLP
[Y,X] = meshgrid(-1:2/N:1);
for tj = 1:length(tvar);
    nl = round(sparsity*(N+1)^2);
    t = 0:1/nl:1;
    t = t(:);
    omega1 = 40*pi;
    omega2 = 35*pi;
    wave_numberx = N*cos(omega1*t).*cos(omega2*t)/2;
    wave_numbery = N*cos(omega1*t).*sin(omega2*t)/2;
    plot(wave_numberx,wave_numbery,'.','MarkerSize',14)
    axis square
    kx = round(N*cos(omega1*t).*cos(omega2*t)/2)+N/2+1;
    ky = round(N*cos(omega1*t).*sin(omega2*t)/2)+N/2+1;
    k = [kx ky];
    rand_ind = unique(k,'rows');
    nl = size(rand_ind,1);
    mask = zeros(N+1);
    for j=1:nl
        mask(rand_ind(j,1),rand_ind(j,2)) = 1;
    end
    R = ifftshift(mask);
    
    Mag_Truth = make_shepp_logan_image(N+1).*exp(R2D2(X,Y,tvar(tj))+i*FreqF(X,Y,tvar(tj)));
    Mag_Truth = ifftshift(Mag_Truth);
    res_x=-1+(0:2*N-1)/N;
    h=res_x(2)-res_x(1);    %delta x
    res_xm = repmat(res_x(:),1,2*N);
    res_ym = repmat(res_x(:)',2*N,1);
    wave_number = -N/2:N/2;
    Mag_Res = make_shepp_logan_image(2*N).*exp(R2D2(res_xm,res_ym,tvar(tj))+i*FreqF(res_xm,res_ym,tvar(tj)));
    F = zeros(N+1);
    h_wait = waitbar(0,'Calculate FFT');
    S_N = length(wave_numberx);
    S = zeros(S_N,1);
    for jk = 1:S_N
        waitbar(jk/S_N,h_wait);
        S(jk)=sum(sum(Mag_Res.*exp(-i*(wave_numberx(jk)*pi*res_xm+wave_numbery(jk)*pi*res_ym))*h^2));
    end
    close(h_wait);
    k_m = round([wave_numberx(:) wave_numbery(:)]); 
    [k_nn,indr,indc] = unique(k_m,'rows');
    F_err = zeros(N+1);
    for jinf = 1:size(k_nn,1)
        cindx = k_nn(jinf,1);
        cindy = k_nn(jinf,2);
        use_ind = find(abs(wave_numberx-cindx)<Window_Patch &abs(wave_numbery-cindy)<Window_Patch);
        Fiterp = scatteredInterpolant(wave_numberx(use_ind),wave_numbery(use_ind),S(use_ind),'nearest');
%         us_ind = find(indc==jinf);
%         fc = mean(S(us_ind,:),1);
        cind = k_nn(jinf,:)+N/2+1;
        F(cind(1,1),cind(1,2)) = Fiterp(cindx,cindy);
%         F(cind(1,1),cind(1,2)) = fc;
%         F(cind(1,1),cind(1,2)) = sum(sum(Mag_Res.*exp(-i*(wave_number(cind(1,1))*pi*res_xm+wave_number(cind(1,2))*pi*res_ym))*h^2));
%         F_err(cind(1,1),cind(1,2)) = fc - F(cind(1,1),cind(1,2));
    end
%     F = F+F_err/2;
    
    F=.5^2*ifftshift(F)*(N+1);
    
    
    
    % build the sampling matrix, R
    
    
    
    
    
    % Recover the image
    % [recovered_TV, ce_TV] = SB(R,F, mu, lambda, gamma_tv,n_inner,n_outer,Mag_Truth);
    
    [recovered_TV, ce_TV] = SB_SPA(R,F, mu, lambda, gamma_tv,n_inner,n_outer,1,Mag_Truth);
    recovered_TV = fftshift(real(recovered_TV));
    recovered_Fourier = fftshift(real(ifft2((N+1)*R.*F)));
    Mag_Truth = fftshift(real(Mag_Truth));
    Er_F=abs(recovered_Fourier-Mag_Truth);
    Er_TV=abs(recovered_TV-Mag_Truth);
    
    clim_max = max(Er_F(:));
    clim_max = max(clim_max,max(Er_TV(:)));
    
    
    h=figure;
    imagesc(fftshift(R));
    axis off square
    colormap gray
    
    h=figure;
    imagesc(Mag_Truth);
    axis off square
    colormap gray
    
    h=figure;
    imagesc(recovered_Fourier);
    axis off square
    colormap gray
    
    disp(['Fourier Err: ',num2str(sqrt(sum(sum(Er_F.^2))))])
    disp(['TV Err: ',num2str(sqrt(sum(sum(Er_TV.^2))))])
    
    h=figure;
    imagesc(recovered_TV);
    axis off square
    colormap gray
    marker = 1;
end



