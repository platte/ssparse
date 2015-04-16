
close all
clear all
clc

%% Sampling Parameters %%
N = 64; % The image will be NxN
sub_res = 4;
nl = 3000;     % Machine number of samples per time window.
N_sub = N*sub_res;
Tstar = 4*pi;
tvar = 1:pi:8*pi;
time_window = .5;
omega1 = 107.1416;
omega2 = 23.1416;
gamma_dB = 40;  %  PSNR Noise to be added in Db
do_noise = 1;

%% Uses SL Phantom and Rstar and Frequency Functions defined below %%
R2D2 = @(x,y,t) -t/Tstar;
FreqF = @(x,y,t) .25*(exp(-12*((y-.6).^2+x.^2))+exp(-12*((y+.6).^2+x.^2)))*t;

%% Split Bregman Optimization Parameters %%
n_inner = 10;
n_outer = 200;
mu = 1;
lambda = 100;
gamma_tv = 1;
m_order = 1;

%% Time Series for image reconstrction %%

Mag_Truth_ts = zeros(N+1,N+1,length(tvar));
Mag_TV_ts = zeros(N+1,N+1,length(tvar));
Mag_Fourier_ts = zeros(N+1,N+1,length(tvar));
h_fig = figure;
for tj = 1:length(tvar);
    %% Spiral Sampling %%
    t = 0:1/nl:1;
    t = t(:);
    wave_numberx = N*cos(omega1*t).*cos(omega2*t)/2;
    wave_numbery = N*cos(omega1*t).*sin(omega2*t)/2;
    figure(h_fig)
    subplot(2,3,1)
    plot(wave_numberx,wave_numbery,'.','MarkerSize',14)
    axis square
    title('Measurement Sampling')
    kx = round(N_sub*cos(omega1*t).*cos(omega2*t)/2)+N_sub/2+1;
    ky = round(N_sub*cos(omega1*t).*sin(omega2*t)/2)+N_sub/2+1;
    k = [kx ky];
    rand_ind = unique(k,'rows');
    nl = size(rand_ind,1);
    mask = zeros(N_sub+1);
    for j=1:nl
        mask(rand_ind(j,1),rand_ind(j,2)) = 1;
    end
    R = ifftshift(mask);
    [Y,X] = meshgrid(-1:2/N:1);
    Mag_Truth_sub = make_shepp_logan_image(N+1).*exp(R2D2(X,Y,tvar(tj))+i*FreqF(X,Y,tvar(tj)));
    Mag_Truth = zeros(N_sub+1);
    Mag_Truth(N_sub/2+1-N/2:N_sub/2+1+N/2,N_sub/2+1-N/2:N_sub/2+1+N/2) = Mag_Truth_sub;
    Mag_Truth = ifftshift(Mag_Truth);
    res_x=-1+(0:2*N_sub-1)/N_sub;
    h=res_x(2)-res_x(1);    %delta x
    res_xm = repmat(res_x(:),1,2*N_sub);
    res_ym = repmat(res_x(:)',2*N_sub,1);
    
    S_N = length(wave_numberx);
    S = zeros(S_N,1);
    h_wait = waitbar(0,'Calculate FFT');
    kx = round(N_sub*cos(omega1*t).*cos(omega2*t)/2);
    ky = round(N_sub*cos(omega1*t).*sin(omega2*t)/2);
    k = [kx ky];
    k = unique(k,'rows')./sub_res;
    
    for jk = 1:S_N
        waitbar(jk/S_N,h_wait);
        Mag_Res = make_shepp_logan_image(2*N_sub).*exp(R2D2(res_xm,res_ym,tvar(tj)+time_window*(2*t(jk)-1))+i*FreqF(res_xm,res_ym,tvar(tj)+time_window*(2*t(jk)-1)));
        S(jk)=sum(sum(Mag_Res.*exp(-i*(wave_numberx(jk)*pi*res_xm+wave_numbery(jk)*pi*res_ym))*h^2));
    end
    close(h_wait);
    %% Add Noise %%
    noise_added_r = randn(size(S));
    noise_added_i = randn(size(S));
    gamma_r = 10^(log10(sqrt(mean(real(S(:)).^2))/sqrt(mean(noise_added_r(:).^2)))-gamma_dB/20);
    gamma_i = 10^(log10(sqrt(mean(imag(S(:)).^2))/sqrt(mean(noise_added_i(:).^2)))-gamma_dB/20);
    S = S+do_noise*(gamma_r*noise_added_r+gamma_i*noise_added_i);
    
    
    Fiterp = scatteredInterpolant(wave_numberx,wave_numbery,S,'nearest');
    Fv = Fiterp(k(:,1),k(:,2));
    F = zeros(N_sub+1);
    for jk = 1:size(k,1)
        F(k(jk,1)*sub_res +N_sub/2+1,k(jk,2)*sub_res +N_sub/2+1) = Fv(jk);
    end
    F=.5^2*ifftshift(F)*(N+1)/sub_res;

    % Recover the image
    % [recovered_TV, ce_TV] = SB(R,F, mu, lambda, gamma_tv,n_inner,n_outer,Mag_Truth);
    
    [recovered_TV, ce_TV] = SB_SPA(R,F, mu, lambda, gamma_tv,n_inner,n_outer,m_order,Mag_Truth);
    recovered_TV = fftshift(recovered_TV);
    recovered_Fourier = fftshift(ifft2((N+1)*R.*F));
    Mag_Truth = fftshift(Mag_Truth);
    Mag_Truth = Mag_Truth(N_sub/2+1-N/2:N_sub/2+1+N/2,N_sub/2+1-N/2:N_sub/2+1+N/2);
    recovered_Fourier = recovered_Fourier(N_sub/2+1-N/2:N_sub/2+1+N/2,N_sub/2+1-N/2:N_sub/2+1+N/2);
    recovered_TV = recovered_TV(N_sub/2+1-N/2:N_sub/2+1+N/2,N_sub/2+1-N/2:N_sub/2+1+N/2);
    
    Mag_Truth_ts(:,:,tj) = Mag_Truth;
    Mag_TV_ts(:,:,tj) = recovered_TV;
    Mag_Fourier_ts(:,:,tj) = recovered_Fourier;
    
    Er_F=abs(recovered_Fourier-Mag_Truth);
    Er_TV=abs(recovered_TV-Mag_Truth);
    
    figure(h_fig)
    subplot(2,3,2)
    imagesc(fftshift(R));
    axis off square
    colormap gray
    title('Optimization Sampling')
    
    figure(h_fig)
    subplot(2,3,3)
    imagesc(real(Mag_Truth));
    axis off square
    colormap gray
    title('Ground Truth')
    
    figure(h_fig)
    subplot(2,3,4)
    imagesc(real(recovered_Fourier));
    axis off square
    colormap gray
    title('Fourier Reconstruction')
    
    disp(['Fourier Err: ',num2str(sqrt(sum(sum(Er_F.^2))./sum(abs(Mag_Truth(:)).^2)))])
    disp(['TV Err: ',num2str(sqrt(sum(sum(Er_TV.^2))./sum(abs(Mag_Truth(:)).^2)))])
    
    figure(h_fig)
    subplot(2,3,5)
    imagesc(real(recovered_TV));
    axis off square
    colormap gray
    title('TV Reconstruction')
    drawnow
end

do_ind = 1:4;
R2D2_TV = zeros(N+1);
FreqF_TV = zeros(N+1);
M_TV = zeros(N+1);
R2D2_Fourier = zeros(N+1);
FreqF_Fourier = zeros(N+1);
M_Fourier = zeros(N+1);
dt = tvar(do_ind);
A = [dt(:).^0 dt(:)];
AT = A';
ATA_inv = inv(AT*A);
for jx = 1:N+1
    for jy = 1:N+1
        ts = Mag_TV_ts(jx,jy,do_ind);
        ts = log(ts(:));
        ls_co = ATA_inv*AT*ts;
        M_TV(jx,jy) = real(exp(ls_co(1)));
        R2D2_TV(jx,jy) = -real(ls_co(2));
        FreqF_TV(jx,jy) = imag(ls_co(2));
        
        ts = Mag_Fourier_ts(jx,jy,do_ind);
        ts = log(ts(:));
        ls_co = ATA_inv*AT*ts;
        M_Fourier(jx,jy) = real(exp(ls_co(1)));
        R2D2_Fourier(jx,jy) = -real(ls_co(2));
        FreqF_Fourier(jx,jy) = imag(ls_co(2));
    end
end
figure(h_fig)
subplot(2,3,1)
imagesc(M_Fourier);
axis off square
colormap gray
colorbar
title('m(x,y) Fourier')
subplot(2,3,2)
imagesc(M_TV);
axis off square
colormap gray
colorbar
title('m(x,y) TV')
subplot(2,3,3)
imagesc(make_shepp_logan_image(N+1));
axis off square
colormap gray
colorbar
title('m(x,y) Ground Truth')

Mask_Fourier(find(abs(Mask_Fourier)<.1)) = 0;
Mask_Fourier(find(abs(Mask_Fourier)>0)) = 1;
subplot(2,3,4)
imagesc(FreqF_Fourier.*Mask_Fourier);
axis off square
colormap gray
colorbar
title('f(x,y) Fourier')

Mask_TV(find(abs(Mask_TV)<.2)) = 0;
Mask_TV(find(abs(Mask_TV)>0)) = 1;
subplot(2,3,5)
imagesc(FreqF_TV.*Mask_TV);
axis off square
colormap gray
colorbar
title('f(x,y) TV')
Mask_GT = make_shepp_logan_image(N+1);
Mask_GT(find(abs(Mask_GT)<.1)) = 0;
Mask_GT(find(abs(Mask_GT)>0)) = 1;
subplot(2,3,6)
imagesc(Mask_GT.*FreqF(X,Y,1));
axis off square
colormap gray
colorbar
title('m(x,y) Ground Truth')

disp(['Tstar: Fourier = ',num2str(1/mean(R2D2_Fourier(:))),' TV = ',num2str(1/mean(R2D2_TV(:))),' Ground Truth= ',num2str(Tstar)])
disp(['Fourier Err f = ',num2str(sqrt(sum(sum((abs(FreqF(X,Y,1)-FreqF_Fourier).*Mask_Fourier).^2))/sum(sum((abs(FreqF(X,Y,1)).*Mask_TV).^2))))]);
disp(['TV Err f = ',num2str(sqrt(sum(sum((abs(FreqF(X,Y,1)-FreqF_TV).*Mask_TV).^2))/sum(sum((abs(FreqF(X,Y,1)).*Mask_TV).^2))))]);

disp(['Fourier Err m = ',num2str(sqrt(sum(sum((abs(make_shepp_logan_image(N+1)-M_Fourier).*Mask_Fourier).^2))/sum(sum((abs(make_shepp_logan_image(N+1)).*Mask_TV).^2))))]);
disp(['TV Err m = ',num2str(sqrt(sum(sum((abs(make_shepp_logan_image(N+1)-M_TV).*Mask_TV).^2))/sum(sum((abs(make_shepp_logan_image(N+1)).*Mask_TV).^2))))]);

