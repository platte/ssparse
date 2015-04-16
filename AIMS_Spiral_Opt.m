
close all
clear all
clc

%% Try to find best sprial %%

N = 64; % The image will be NxN
nl = 3000;
sub_res = 4;


N_sub = N*sub_res;
omega1 = pi:2:100*pi;
omega2 = pi:2:100*pi;


t = 0:1/nl:1;
t = t(:);
num_points = zeros(length(omega1),length(omega2));
round_off_err = zeros(length(omega1),length(omega2));
h_wait = waitbar(0,'Testing Spirals');
for jo1=1:length(omega1)
    waitbar(jo1/length(omega1),h_wait);
    for jo2=1:length(omega2)
        %% Create measured and interpolated sample points %%
        wave_numberx = N_sub*cos(omega1(jo1)*t).*cos(omega2(jo2)*t)/2;
        wave_numbery = N_sub*cos(omega1(jo1)*t).*sin(omega2(jo2)*t)/2;
        
        kx = round(wave_numberx);
        ky = round(wave_numbery);
        k = [kx ky];
        rand_ind = unique(k,'rows');
        %% Determine number of points used in reconstruction %%
        num_points(jo1,jo2) = size(rand_ind,1);
        
        %% Estimate interpolation error %%
        Fiterp = scatteredInterpolant(wave_numberx,wave_numbery,wave_numberx+wave_numbery,'nearest');
        Fv = Fiterp(k(:,1),k(:,2));
        FvT = k(:,1)+k(:,2);
        round_off_err(jo1,jo2) = sum(abs(Fv-FvT));
    end
end
close(h_wait)
h_fig = figure;
figure(h_fig)
subplot(1,3,1)
imagesc(num_points)
axis off square
colormap gray
title('Number of Unique Points')
subplot(1,3,2)
imagesc(round_off_err)
axis off square
colormap gray
title('Round off error')
subplot(1,3,3)
imagesc(num_points./round_off_err)
axis off square
colormap gray
title('Round off error')

comb_measure = num_points./round_off_err;
max_v = max(comb_measure(:));

[max_ind1,max_ind2] = find(comb_measure==max_v);
disp(['Best \omega_1: ',num2str(omega1(max_ind1)),' Best \omega_2: ',num2str(omega2(max_ind2))])

