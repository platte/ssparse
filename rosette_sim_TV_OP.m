%rosette_sim.m
%Edward Walsh, Ph.D.
%OBJECTIVE: Computes a variable density rosette trajectory with both slew rate
%           and amplitude limiting.  The angular frequencies of oscillation
%           and acquisition duration can be specified.  Output signal is in
%           array named 'sig'. Program loads file 'phantom1.mat' which contains
%           array 'E' which specified the Shepp-Logan phantom image.  The
%           first six columns are the regular parameter array for the
%           Shepp-Logan phantom ellipse specifications.  Column 7 are the
%           R2* values for each ellipse, and Column 8 is the frequency
%           offset for the elipses.
%VARIABLE DEFINITIONS:  FOV=field of view
%                       Nres=image resolution
%                       Nres_samp=size of sampling grid
%                       N_shots=number of interleaves
%                       phi=angular step between interleaves
%                       pt1,pt2=array indices for final image in large FOV
%                       E=array of ellipse specifications,R2*,frequency
%                       n_points=number of sample points if fully sampled
%                       om1=angular frequency of oscillation through k=0
%                       om2=angular frequency of rotation around k=0                       
%                       gamma=gyromagnetic ratio (Hz/G);
%                       point_trans=point # for change to undersampling
%                       kmax=maximum spatial frequency
%                       cent=center coordinate of sampling grid
%                       k_step_grid=spatial frequency step on sample grid
%                       k_step=spatial frequency step per turn
%                       k=basic k-space trajectory
%                       k2=trajectory for interleaves
%                       n,nn,m,x,y,count=loop counters
%                       temp1,temp2=temporary holding variables
%                       sz=signal and k-vector size
%                       Gmax=gradient amplitude limit (G/cm)
%                       Smax=gradient slew rate limit (G/cm/s)
%                       Tacq=acquisition time for limited trajectory
%                       nt=time vector
%                       g=gradient signal for limited (at end) trajectory
%                       s=slew signal for limited (at end) trajectory
%                       Gmax_lim,Smax_lim=max values in limited signals
%                       sig=simulated MR signal for phantom image
%                       sig_current=signals for individual ellipses
%                       sig_image=interleave signals for head phantom
%                       sig_psf=interleave signals for psf
%                       inputsize=size of E
%                       n_ellipse=number of ellipses in phantom image
%                       ellipses=specification array for 'phantom' function
%                       R2_star,offres=vectors of R2*, frequency
%                       I=ellipse intensity
%                       A,B=ellipse axes
%                       x_off,y_off=ellipse position offsets
%                       alpha=ellipse rotation angle
%                       u,v=rotation matrix for signal computation
%                       temp1,temp2=temporary holding variables
%                       decay=R2* decay component of signal
%                       f_offset=off-resonance component of signal         

close all
clear all
clc

%% Hard Code input %%
FOV= 19.2;% input('Field of View (cm) ? =');
Nres=64; %input('Resolution ? =');
Tacq=.05; %input('Acquisition Time ? =');
om1=3576.7; %input('Omega 1 ? =');
om2=3967.9; %input('Omega 2 ? =');
t_dwell=5e-6; %input('Dwell time per point ? =');
load phantom1

%% Simulate Signal %%
N_shots=1;
kmax=Nres./(2.*FOV);
nt=(0:t_dwell:Tacq)';
n_points=max(size(nt));
gamma=4258;
k_lim=kmax.*sin(om1.*nt).*exp(1i.*om2.*nt);
phi=(2.*pi.*(om2./om1))./N_shots;
sz=n_points;
figure(1)
plot(real(k_lim),imag(k_lim));
xlabel('kx (1/cm)');
ylabel('ky (1/cm)');
title('Interleave 1, Amplitude/Slew Limited Trajectory');
g=diff(k_lim)./(gamma.*t_dwell);
figure(2)
plot(nt(1:(n_points-1)),real(g)./100);
xlabel('Time (s)');
ylabel('Amplitude (T/m)');
title('Interleave 1, X gradient component');
s=diff(g)./t_dwell;
figure(3)
plot(nt(1:(n_points-2)),real(s)./100);
xlabel('Time (s)');
ylabel('X slew rate (T/m/s)');
title('Interleave 1, X slew rate');
Gmax_lim=max(real(g));
Smax_lim=max(real(s));
fprintf('Maximum amplitude of limited trajectory = %f T/m.\n',Gmax_lim./100);
fprintf('Maximum slew rate of limited trajectory = %f T/m/s.\n',Smax_lim./100);
fprintf('Readout duration = %f s.\n',max(nt));

%ASSEMBLE SIGNALS FROM ELLIPSES
inputsize=size(E);
n_ellipse=inputsize(1);
R2_star=E(1:n_ellipse,7);
offres=E(1:n_ellipse,8);
sig_image(1:sz,1:N_shots)=0;
k2 = zeros(sz,N_shots);
k1 = zeros(sz,N_shots);
kx1 = zeros(sz,N_shots);
ky1 = zeros(sz,N_shots);
for m=1:N_shots
    sig = zeros(sz,1);
    sig_current(1:sz,1:n_ellipse)=0+0.*1i;
    %k-components
    sgn=-1.^m;
    k2(:,m)=k_lim.*exp(sgn.*1i.*(m-1).*phi);
    kx=real(k2(:,m));
    ky=imag(k2(:,m));
    kx1(:,m)=kx;
    ky1(:,m)=ky;
    for n=1:n_ellipse
        %extract ellipse parameters and convert normalized values to FOV
        %x,y directions and signs to permit use of unmodified 'phantom' arrays
        I=E(n,1);
        B=E(n,2).*FOV./2;
        A=E(n,3).*FOV./2;
        y_off=-E(n,4).*FOV./2;
        x_off=E(n,5).*FOV./2;
        alpha=E(n,6).*(pi/180);
        %frequency space representation of ellipse
        u=kx.*cos(alpha)+ky.*sin(alpha);
        v=-kx.*sin(alpha)+ky.*cos(alpha);
        temp1=sqrt((u.*A./B).^2+v.^2);
        temp2=2.*pi.*B.*temp1;
        sig_current(:,n)=I.*exp(1i.*2.*pi.*(kx.*x_off+ky.*y_off)).*(A.*besselj(1,temp2)./(temp1+0.00001));
        decay=exp(-R2_star(n).*nt);
        f_offset=exp(1i.*2.*pi.*offres(n).*nt);
        sig_current(:,n)=sig_current(:,n).*decay.*f_offset;
        sig=sig+sig_current(:,n);
    end%loop over ellipses
    sig_image(:,m)=sig;
end%loop over shots
%FIDIN=zeros(sz-1,2,2,2);%for compatibility with rundtb3sim.m
%FIDIN(1:sz-1,1,1,1)=sig(1:sz-1);
figure(4)
plot(nt,abs(sig));
title('Signal Magnitude')
xlabel('Tacq (s)')
kss=k_lim(1:sz-2);
clear phi
phi(1:sz-2)=0;%residual current compensation not relevant for simulations
%save Ksim kss sig

%% Nearest Neighbor Interpolation to Uniform K-Space %%
Nk = 12*kmax*FOV/2; %round(sqrt(n_points))/2;
scaling_factor = Nk/kmax;

Fiterp = scatteredInterpolant(scaling_factor*kx,scaling_factor*ky,sig_image,'nearest');
k = unique(round([scaling_factor*kx(:) scaling_factor*ky(:)]),'rows');
Fv = Fiterp(k(:,1),k(:,2));
F = zeros(2*Nk+1);
R = zeros(2*Nk+1);


for jk = 1:size(k,1)
    F(k(jk,1)+Nk+1,k(jk,2)+Nk+1) = Fv(jk);
    R(k(jk,1)+Nk+1,k(jk,2)+Nk+1) = 1;
end

F(end,:) = [];
F(:,end) = [];
R(end,:) = [];
R(:,end) = [];

sub_res = 1/(kmax*FOV/Nk);
Fp = zeros(2*sub_res*Nk);
Rp = zeros(2*sub_res*Nk);
Fp(sub_res*Nk-Nk+1:sub_res*Nk+Nk,sub_res*Nk-Nk+1:sub_res*Nk+Nk) = F;
Rp(sub_res*Nk-Nk+1:sub_res*Nk+Nk,sub_res*Nk-Nk+1:sub_res*Nk+Nk) = R;

Image_rec = (2*sub_res*Nk)^2*ifftshift(ifft2(fftshift(Fp)));
Image_rec = real(Image_rec(sub_res*Nk-Nk+1:sub_res*Nk+Nk,sub_res*Nk-Nk+1:sub_res*Nk+Nk));


%% TV optimization %%

%% Split Bregman Optimization Parameters %%
n_inner = 5;
n_outer = 100;
mu = 1;
lambda = 100;
gamma_tv = 0;
m_order = 1;

F = fftshift(Fp);
R = fftshift(Rp);
Image_rec_TV = SB_SPA_MD(R,F, mu, lambda, gamma_tv,n_inner,n_outer,m_order);
Image_rec_TV = (2*sub_res*Nk)^2*ifftshift(Image_rec_TV);
Image_rec_TV = real(Image_rec_TV(sub_res*Nk-Nk+1:sub_res*Nk+Nk,sub_res*Nk-Nk+1:sub_res*Nk+Nk));
figure,
subplot(1,3,1)
imagesc(Image_rec)
axis square off
colormap gray
title('Fourier Reconstruction','FontSize',16)
subplot(1,3,2)
imagesc(Rp(sub_res*Nk-Nk+1:sub_res*Nk+Nk,sub_res*Nk-Nk+1:sub_res*Nk+Nk))
axis square off
colormap gray
title('Nearest Neighbor Reconstruction','FontSize',16)
subplot(1,3,3)
imagesc(Image_rec_TV)
axis square off
colormap gray
title('TV Reconstruction','FontSize',16)




















%%%%%%%% Temp Code to Figure out FFT padding %%%%%%%%%%%%%%%

% [Y,X]=meshgrid(-FOV:FOV/Nk:FOV);
% 
% Forg = zeros(size(X));
% hig_w = waitbar(0,'Calculating');
% for jx = 1:size(X,1);
%     waitbar(jx/size(X,1))
%     for jy=1:size(X,2);
%         Forg(jx,jy) = real(sum(sig_image.*exp(i*(kx*pi*X(jx,jy)+ky*pi*Y(jx,jy)))));
%     end
% end
% close(hig_w);
% figure, 
% subplot(1,3,1)
% imagesc(Fint)
% axis square off
% colorbar
% subplot(1,3,2)
% imagesc(FU)
% axis square off
% colorbar
% subplot(1,3,3)
% imagesc(Forg)
% axis square off
% colorbar
% drawnow
% %% How to reconstruct fft on -1<=x<=1-1/N %%
% N = 4;
% fx= rand(N,1);
% wn = exp(-2*pi*i/N);
% fftx = zeros(N,1);
% for jk=1:N
%     for j=1:N
%         fftx(jk) = fftx(jk)+fx(j)*wn^((j-1-N/2)*(jk-1-N/2)/2);
%     end
% end
% [fft(fx) fftshift(fftx).*(wn.^((-N/2)*((0:N-1)'-N/2)))]

% %% How to reconstruct fft on sub resolution in K on 0=x<=1-1/N%%
% sub_res = 3;
% N = 4;
% fx = zeros(sub_res*N,1);
% fx(1:N) = rand(N,1);
% wn = exp(-2*pi*i/N);
% fftx = zeros(N,1);
% 
% for jk=1:N
%     for j=1:N
%         fftx(jk) = fftx(jk)+fx(j)*wn^((j-1)*(jk-1)/sub_res);
%     end
% end
% fsb = fft(fx);
% [fsb(1:N) fftx]

% %% How to reconstruct fft on sub resolution in K on -1<=x<=1-1/N %%
% sub_res = 4;
% N = 6;
% kxu = -N/2:N/2-1;
% XU =(-N/2:N/2-1)/sub_res;
% fx = zeros(sub_res*N,1);
% fx(1:N) = rand(N,1);
% wn = exp(-2*pi*i/N);
% fftx = zeros(N,1);
% 
% for jk=1:N
%     for j=1:N
%         fftx(jk) = fftx(jk)+fx(j)*wn^(XU(j)*kxu(jk));
%     end
% end
% fsb = fftshift(fft(fx));
% [fsb(sub_res*N/2-N/2+1:sub_res*N/2+N/2) fftx.*wn.^((N/2)*((0:N-1)'-N/2)/sub_res)]


%% How to reconstruct fft on sub resolution in K on -1<=x<=1-1/N 2D%%

% sub_res = 2;
% N = 2;
% [kyu,kxu] = meshgrid(-N:N-1);
% [YU,XU]=meshgrid((-N:N-1)/(2*sub_res));
% fx = zeros(2*sub_res*N);
% fx(1:2*N,1:2*N) = rand(2*N,2*N);
% wn = exp(-2*pi*i/N);
% fftx = zeros(2*N);
% 
% for jkx=1:2*N
%     for jky=1:2*N
%         for jx=1:2*N
%             for jy=1:2*N
%                 fftx(jkx,jky) = fftx(jkx,jky)+fx(jx,jy)*wn^(XU(jx,jy)*kxu(jkx,jky)+YU(jx,jy)*kyu(jkx,jky));
%             end
%         end
%     end
% end
% fftx= fftx;
% fsb = fftshift(fft2(fx));
% fsb = fsb(sub_res*N-N+1:sub_res*N+N,sub_res*N-N+1:sub_res*N+N).*(wn.^(-N*kxu/(2*sub_res)-N*kyu/(2*sub_res)));
% [fsb-fftx]

%% How to reconstruct ifft on sub resolution in K on -1<=x<=1-1/N %%

% sub_res = 2;
% N = 2;
% fk = zeros(2*sub_res*N,1);
% fk(1:2*N) = rand(2*N,1);
% kxu = (-N:N-1)';
% xu = (-N:N-1)'/(2*sub_res);
% wn = exp(-2*pi*i/N);
% sifft_fx = zeros(2*N,1);
% for jx=1:2*N
%     for jk=1:2*N
%         sifft_fx(jx) = sifft_fx(jx)+fk(jk)*wn^(-xu(jx)*kxu(jk))/(2*N);
%     end
% end
% 
% ifft_fx = ifftshift(ifft(fk));
% ifft_fx = sub_res*ifft_fx(sub_res*N-N+1:sub_res*N+N,1).*(wn.^(N*kxu/(2*sub_res)));
% ifft_fx-sifft_fx
% 
% %% use zero padding and shift %%
% fk = fk(1:2*N);
% fkp = zeros(2*sub_res*N,1);
% fkp(sub_res*N-N+1:sub_res*N+N,1) = fk;
% 
% ifft_fx = sub_res*ifftshift(ifft(fftshift(fkp)));
% ifft_fx = ifft_fx(sub_res*N-N+1:sub_res*N+N,1);
% ifft_fx-sifft_fx

%% How to reconstruct ifft on sub resolution in K on -1<=x<=1-1/N 2D%%

% sub_res = 2;
% N = 4;
% [kyu,kxu] = meshgrid(-N:N-1);
% [YU,XU]=meshgrid((-N:N-1)/(2*sub_res));
% fx = rand(2*N,2*N);
% wn = exp(-2*pi*i/N);
% fftx = zeros(2*N);
% 
% 
% for jx=1:2*N
%     for jy=1:2*N
%         fftx(jx,jy) = sum(sum(fx.*(wn.^(-XU(jx,jy)*kxu-YU(jx,jy)*kyu))));
%     end
% end
% fxp = zeros(2*sub_res*N);
% fxp(sub_res*N-N+1:sub_res*N+N,sub_res*N-N+1:sub_res*N+N) = fx;
% 
% fsb = ifftshift(ifft2(fftshift(fxp)));
% fsb = (2*sub_res*N)^2*fsb(sub_res*N-N+1:sub_res*N+N,sub_res*N-N+1:sub_res*N+N);
% [fsb-fftx]




% sub_res = 2;
% N = 2;
% [kyu,kxu] = meshgrid(-N:N-1);
% [YU,XU]=meshgrid((-N:N-1)/(2*sub_res));
% fx = zeros(2*sub_res*N);
% fx(1:2*N,1:2*N) = rand(2*N,2*N);
% wn = exp(2*pi*i/N);
% fftx = zeros(2*N);
% 
% for jkx=1:2*N
%     for jky=1:2*N
%         for jx=1:2*N
%             for jy=1:2*N
%                 fftx(jkx,jky) = fftx(jkx,jky)+fx(jx,jy)*wn^(XU(jx,jy)*kxu(jkx,jky)+YU(jx,jy)*kyu(jkx,jky));
%             end
%         end
%     end
% end
% fftx= fftx;
% fsb = ifftshift(ifft2(fx));
% fsb = (2*sub_res*N)^2*fsb(sub_res*N-N+1:sub_res*N+N,sub_res*N-N+1:sub_res*N+N).*(wn.^(-N*kxu/(2*sub_res)-N*kyu/(2*sub_res)));
% [fsb-fftx]

