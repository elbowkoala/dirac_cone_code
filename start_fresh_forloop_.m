tic;
%name = '20170518_00035.fits';
%info = fitsinfo(name);

load 'purply_colormap_64.mat';
load 'normed_binned_cones.mat';

%% Set Parameters 
pix2eV = (1.599/(2*496));
pix2invA = 0.512*0.04631/180*3.1415*14/30*sqrt(110-4);

X_size = 31;
Y_size = 31;
num_scans = X_size *Y_size;
X = zeros(1, num_scans);
Y = zeros(1, num_scans);

frame_size = [768,997];
cone_range_K=[351,700];
cone_range_E=[1,800];
e_bin = 3;
k_bin = 1;

%DPI_K_range = round(-50/k_bin):round(50/k_bin);
DPI_E_range = round(300/e_bin):round(450/e_bin);

first_run=0;
cluster_radius=1;

%% First Run Resets All Maps & Calibrates CCD Pixels by Column(Energy) 
%%  & Stores the Blob-Normalized, Binned, (& Clustered, if cluster_radius>0) Cones 
if first_run==1    
    load('P pol data.mat');
       
    dirac_new_map = zeros(1,num_scans);
    dirac_E_map = zeros(1,num_scans);
    dirac_K_map = zeros(1,num_scans);
    dirac_V_map = zeros(1,num_scans);
    fit_evals_map = zeros(1,num_scans);
    DP_intensity_map = zeros(1,num_scans);

    blob_row_range = (1:200); %%frame pixel row range to inspect sparse blobs 
    [blob_I_normed_cones, ave_blob_I_x, ave_blob_I_y] ...
        = blob_I_equalizer(data,frame_size,num_scans,cone_range_K,cone_range_E,blob_row_range,50,1,850);
    
    get_binned_size = Binning_2d(blob_I_normed_cones(:,:,1),e_bin,k_bin);
    binned_cones = zeros(size(get_binned_size,1),size(get_binned_size,2),num_scans);
    raw_full_binned_cone = zeros(size(get_binned_size));
    if cluster_radius > 0
        cluster_binned_cones = zeros(size(get_binned_size,1),size(get_binned_size,2),num_scans);
    end
    for i = 1:num_scans
        if rem(i,50)==0
            disp(i)
        end
        X(i) = data{1,2}(i);
        Y(i) = data{1,3}(i);
        
        cone = blob_I_normed_cones(:,:,i);
        binned_cone = Binning_2d(cone,e_bin,k_bin);       
        binned_cones(:,:,i) = binned_cone;
        raw_full_binned_cone = raw_full_binned_cone + binned_cone;
        DP_intensity_map(i) = sum(sum(binned_cone(:,DPI_E_range)));
        
        if cluster_radius > 0     
            [i_row,i_col] = ind2sub([31,31],i);
            to_add = zeros(31);
            %cone = zeros(length(cone_range_K(1):cone_range_K(end)),length(cone_range_E(1):cone_range_E(end)));
            for i_r = i_row-2 : i_row+2
                if i_r<1 || i_r>31
                    continue
                end
                for i_c = i_col-2 : i_col+2
                    if i_c<1 || i_c>31
                        continue
                    end
                    dist_from_i = sqrt((i_r-i_row)^2 + (i_c-i_col)^2);
                    if dist_from_i > 0 && dist_from_i <= cluster_radius
                        to_add(i_r,i_c) = 1;
                    end
                end
            end
            for jjj= find(to_add==1)'
                cone = cone + blob_I_normed_cones(:,:,jjj);
            end
            cluster_num = length(find(to_add==1));
            cone = cone / (1+cluster_num);
    
            bcone = Binning_2d(cone,e_bin,k_bin);
            cluster_binned_cones(:,:,i) = bcone;
        end     
    end
    DPI_map = reshape(DP_intensity_map,[X_size,Y_size]);
    X=reshape(X,[X_size,Y_size]);
    Y=reshape(Y,[X_size,Y_size]);

    figure, imagesc(Y(1,:),X(:,1),DPI_map); axis xy; title('DP Intensity Map');
    xlabel('Y (mm)');
    ylabel('X (mm)');   
end

DPI_mean = mean(DP_intensity_map(:));
DPI_std = std(DP_intensity_map(:));

%%
kLOS_Ilow_thresh = 0;
kLOS_Iup_thresh = 4;
%draw_y_range = round(250/e_bin):1:round(550/e_bin);
trial_DEs = round(280/e_bin):1:round(450/e_bin);
trial_DVs = (1/e_bin)*(pix2invA/pix2eV)* (2.5:0.1:5.0);  %Expect a slope of ~3.5eV-Ang
draw_kgamma = 4;  draw_egamma = 5;
draw_bel_DE = 50; draw_abo_DE = 50;
fade_bel_DE = 15; fade_abo_DE = 35;
leg_fade = zeros(1,size(binned_cones(:,:,1),2));
SZ2 = round(size(binned_cones(:,:,1),2)/2);
for i = SZ2-draw_bel_DE:SZ2 - fade_bel_DE
leg_fade(i) = (i-(SZ2-draw_bel_DE))/(draw_bel_DE-fade_bel_DE);
end
for i = SZ2+fade_abo_DE:SZ2+draw_abo_DE
leg_fade(i) = (i-(SZ2+draw_abo_DE))/(fade_abo_DE-draw_abo_DE);
end
leg_fade(SZ2-fade_bel_DE:SZ2+fade_abo_DE) = 1;
leg_fader = repmat(leg_fade,size(binned_cones(:,:,1),1),1);

%init_kLOS_map = zeros(1,961);
plot_fig=0;
do_rfc_map=0;
draw_DV_matrices = XVFg_matrices_drawer2(trial_DVs,draw_kgamma,draw_egamma,size(binned_cones(:,:,1)));
tic
for i = 813:961%[round(961*rand),round(961*rand),round(961*rand),round(961*rand),round(961*rand)]%1:1:num_scans
    disp(['On scan ',num2str(i)])
    if rem(i,50)==0
        toc
    end
    
    bcone = binned_cones(:,:,i);    
    
    if DP_intensity_map(i) <  DPI_mean-DPI_std
        kLOS_erange = round((1/e_bin)*(100:700));
        kLOS_sigma = 8;
    else
        kLOS_erange = round((1/e_bin)*(200:550));
        kLOS_sigma = 4;
    end
    kLOS_pt = the_kLOS_finder( bcone, kLOS_sigma, kLOS_erange, kLOS_Ilow_thresh, kLOS_Iup_thresh);
    init_kLOS_map(i) = kLOS_pt;         
%     ksym_range = round((1/k_bin)*(-100:100));
%     kaxis_pre_sym = kLOS_pt + ksym_range;    
%     [sym_bcone, sym_kaxis] = Symmetrized_spectra(bcone(kaxis_pre_sym,:),ksym_range');
%     corr_cone = sym_bcone;
%     trial_DKs = round(size(corr_cone,1)/2) + (-5:5);
   
    corr_cone = bcone;

    trial_DKs = kLOS_pt + (-5:5);   
     
    I__of_DVs_table = cell(length(trial_DKs),length(trial_DEs));
    best__DVs_table = NaN*ones(length(trial_DKs),length(trial_DEs));
    max__corrs_table = NaN*ones(length(trial_DKs),length(trial_DEs));
    for trial_DE_i = 1:size(trial_DEs,2)
        trial_DE = trial_DEs(trial_DE_i); 
        
        draw_y_range = trial_DE-draw_bel_DE:trial_DE+draw_abo_DE;
        corr_cone_crop = corr_cone(:,draw_y_range);
        ccc_mean = mean2(corr_cone);
        ccc_denom = sqrt(sum(sum((corr_cone_crop-ccc_mean).^2)));
        
        for trial_DK_i = 1:length(trial_DKs)
            if sqrt( (trial_DK_i-length(trial_DKs)/2)^2 + (trial_DE_i - length(trial_DEs)/2)^2) > 20
                continue
            end
            
            trial_DK = trial_DKs(trial_DK_i); 
            I__of_DVs = zeros(size(trial_DVs));
            for trial_DV_i = 1:length(trial_DVs)                   
                [X_matrix] = X_editor3(leg_fader.*draw_DV_matrices(:,:,trial_DV_i),trial_DE,trial_DK,draw_y_range);
                Xm_mean = mean2(X_matrix);%mean(X_matrix(X_matrix>.0001));
                
                %corr_cone_crop = corr_cone(:,draw_y_range);
                %ccc_mean = mean(corr_cone_crop(:));
                NCC_denom = sqrt(sum(sum((X_matrix-Xm_mean).^2))) * ccc_denom;

                I__of_DVs(trial_DV_i) = sum(sum((X_matrix - Xm_mean).*(corr_cone_crop - ccc_mean)))/NCC_denom;
    
            end               
            best__DV_i = find(I__of_DVs==max(I__of_DVs),1,'first');
            I__of_DVs_table{trial_DK_i,trial_DE_i} = I__of_DVs;
            best__DVs_table(trial_DK_i,trial_DE_i) = trial_DVs(best__DV_i);
            max__corrs_table(trial_DK_i,trial_DE_i) = max(I__of_DVs(:));
        end
    end

%     [best_DK_i,best_DE_i] = find(max_corrs_table==nanmax(max_corrs_table(:)));
%     best_I_of_DVs = (I_of_DVs_table{best_DK_i,best_DE_i});
%     best_DK = trial_DKs(best_DK_i);
%     best_DE = trial_DEs(best_DE_i);
%     best_DV_i = find(best_I_of_DVs==max(best_I_of_DVs),1,'first');
%     best_DV = trial_DVs(best_DV_i);

    [best__DK_i,best__DE_i] = find(max__corrs_table==nanmax(max__corrs_table(:)));
    best__I_of_DVs = (I__of_DVs_table{best__DK_i,best__DE_i});
    best__DK = trial_DKs(best__DK_i);
    best__DE = trial_DEs(best__DE_i);
    best__DV_i = find(best__I_of_DVs==max(best__I_of_DVs),1,'first');
    best__DV = trial_DVs(best__DV_i);

    dirac_V_map(i) = best__DV*(e_bin/k_bin);
    dirac_E_map(i) = best__DE*(e_bin);
    dirac_K_map(i) = best__DK*(k_bin);
    fit_evals_map(i) = nanmax(max__corrs_table(:));
    
    if plot_fig == 1
        figure,
        draw_y_range = best__DE-draw_bel_DE:best__DE+draw_abo_DE;
        overlay_matrix = X_editor3(leg_fader.*draw_DV_matrices(:,:,best__DV_i),best__DE,best__DK,draw_y_range);
        overlay_matrix = padarray(overlay_matrix,[0,draw_y_range(1)-1],'pre');
        overlay_matrix = padarray(overlay_matrix,[0,size(corr_cone,2)-draw_y_range(end)],'post');
   
        plotfitx = (1:size(corr_cone,1));
        %figure, subplot(241), imagesc(max_corrs_table), axis xy, title('MC table')
        %subplot(245), imagesc(best_DVs_table), axis xy, title('MC DVs table')
        subplot(242), imagesc(max__corrs_table), axis xy, title('NC MC table')
        subplot(246), imagesc(best__DVs_table), axis xy, title('NC DVs table')
        subplot(2,4,[3,7]), imagesc((1:size(corr_cone,1)),fliplr(1:size(corr_cone,2)),rot90(imgaussfilt(corr_cone,3)+.3*(mean2(imgaussfilt(corr_cone,3))/mean2(overlay_matrix))*overlay_matrix)), axis xy, hold on;

        subplot(244), 
        trial_DVs_units = trial_DVs*(e_bin)*(pix2eV/pix2invA);
        plot(trial_DVs_units,best__I_of_DVs,'b'),  
        xlim([trial_DVs_units(1),trial_DVs_units(end)])
        subplot(248), imagesc((1:size(corr_cone,1)),fliplr(1:size(corr_cone,2)),rot90(imgaussfilt(corr_cone,3))), axis xy, hold on;
        plot(plotfitx, -best__DV*(plotfitx-best__DK)+best__DE, 'c:'), hold on;
        plot(plotfitx, +best__DV*(plotfitx-best__DK)+best__DE, 'c:')

        suptit = suptitle({['i=',num2str(i)];...
            ['NC: DE=',num2str(best__DE),' DV=',num2str(best__DV*(e_bin/k_bin)*(pix2eV/pix2invA)),' MaxVal=',num2str(fit_evals_map(i))]});
        set(suptit,'FontSize',8,'FontWeight','normal')
        pause(.05)
    end

    
end
toc

%% Compare to normxcorr2 method
if do_rfc_map == 1
    
    raw_full_binned_cone = zeros(size(binned_cones(:,:,1)));
    for i = 1:961
        raw_full_binned_cone = raw_full_binned_cone + binned_cones(:,:,i);
    end
    nxc2_Es = zeros(1,961);
    nxc2_Ks = zeros(1,961);
    the_temp = raw_full_binned_cone(100:250,round(260/e_bin):round(500/e_bin));
    rfbc_Kpt = 83; rfbc_Ept = 33;
    for i = 1:961
    the_img = imgaussfilt(binned_cones(:,:,i),5);
    nxc2_scan = normxcorr2(the_temp,the_img);
    [Kpk, Epk] = find(nxc2_scan==max(nxc2_scan(:)));
    K_offset = Kpk-size(the_temp,1);
    E_offset = Epk-size(the_temp,2);
    nxc2_Es(i) = (E_offset + rfbc_Ept + 1)*e_bin;
    nxc2_Ks(i) = (K_offset + rfbc_Kpt + 1)*k_bin;
    end
    %figure, subplot(121), imagesc(reshape(nxc2_Es,31,31)), axis xy, caxis([340,400])
    %subplot(122), imagesc(reshape(nxc2_Ks,31,31)), axis xy
end



%%
figure,
subplot(2,3,1), 
imagesc(imgaussfilt(reshape(pix2eV*(rfc_FL_Es - dirac_E_map),[X_size,Y_size]),.6)); axis xy, set(gca, 'xtick', []), set(gca, 'ytick',[]); 
c1=colorbar('FontSize',8); c1.Label.String='Dirac Energy (abs. eV)';
%caxis(pix2eV*(mean(rfc_FL_Es)-[mean(dirac_E_map)+1.5*std(dirac_E_map), mean(dirac_E_map)-1.5*std(dirac_E_map)]))
caxis([.42,.52])
title('DrawX Method: DP Energy Map');  

subplot(2,3,4)
imagesc(imgaussfilt(reshape((pix2eV/pix2invA)*dirac_V_map,[X_size,Y_size]),.6)); axis xy; set(gca, 'xtick', []), set(gca, 'ytick',[]);
c4=colorbar('FontSize',8); c4.Label.String='Dirac Velocity (eV-A)'; 
caxis([3,4])
title('DrawX Method: DP Vel. Map');

subplot(2,3,2)
imagesc(reshape(pix2invA*(dirac_K_map-mean(dirac_K_map)),[X_size,Y_size])); axis xy; set(gca, 'xtick', []), set(gca, 'ytick',[]);
c2=colorbar('FontSize',8); c2.Label.String='Relative Momentum (inv.A)';
caxis([-.04, .04]);
title('DrawX Method: DP rel. K Map');

subplot(2,3,3)
imagesc(reshape(DP_intensity_map,[X_size,Y_size])); axis xy; set(gca, 'xtick', []), set(gca, 'ytick',[]);
c3=colorbar('Location','eastoutside','FontSize',8); c3.Label.String='DP Region Intensity (A.U.)'; 
title('DrawX Method: DP Region Int. Map');

subplot(2,3,5)
imagesc(reshape(pix2eV*(rfc_FL_Es-nxc2_Es),[X_size,Y_size])); axis xy; set(gca, 'xtick', []), set(gca, 'ytick',[]);
c5=colorbar('FontSize',8); c5.Label.String='Dirac Energy (eV)';
%caxis(pix2eV*(mean(rfc_FL_Es)-[mean(nxc2_Es)+1.5*std(nxc2_Es), mean(nxc2_Es)-1.5*std(nxc2_Es)]))
caxis([.4,.5])
title('RFC Method: DP Energy Map')

subplot(2,3,6)
imagesc(Y(1,:),X(:,1),reshape(pix2invA*(nxc2_Ks-mean(nxc2_Ks)),[X_size,Y_size])); axis xy; %set(gca, 'xtick', []), set(gca, 'ytick',[]);
c6=colorbar('FontSize',8); c6.Label.String='Relative Momentum (inv.A)'; 
caxis([-.04,.04])
title('RFC Method: DP Rel. K')
%imagesc(Y(1,:),X(:,1),dirac_map);

suptitle('Scan 18/02/27: NormCorr - Cluster4')
colormap(map_of_color_64)



% figure,
% subplot(121), imagesc(reshape(DP_intensity_map,31,31)), axis xy
% subplot(122), imagesc(reshape(dirac_K_map,31,31)), axis xy
% 
%{


%}

%     
%     bcone_dos = sum(repmat(abs(-kLOS_pt+[1:size(bcone_mod,1)]'),[1,size(bcone_mod,2)]).*bcone_mod,1);
%     
%     ksym_range = round((1/k_bin)*(-100:100));
%     kaxis_pre_sym = kLOS_pt + ksym_range;    
%     [sym_bcone, sym_kaxis] = Symmetrized_spectra(bcone_mod(kaxis_pre_sym,:),ksym_range');%1:size(cone_mod,1)]');
%     ksym_dos = sum(repmat(abs(sym_kaxis),[1,size(sym_bcone,2)]).*sym_bcone,1);
%        
%    
%         figure, 
%     subplot(221), imagesc(kLOS_pt+[-100:100],fliplr(1:size(bcone_mod,2)),rot90(bcone_mod(kLOS_pt+[-100:100],:))), axis xy, hold on;
%     plot([kLOS_pt,kLOS_pt],[1,size(bcone_mod,2)],'w'), 
%     title(['DPI=',num2str(DP_intensity_map(1,i))])
%     subplot(223), plot(bcone_dos)   
%     
%     subplot(222), imagesc(sym_kaxis,fliplr(1:size(sym_bcone,2)),rot90(sym_bcone)), axis xy, hold on;
%     plot([0,0],[1,size(sym_bcone,2)],'w')
%     title(['DPI=',num2str(DP_intensity_map(2,i))])
%     subplot(224), plot(ksym_dos)
%     suptitle(['i=',num2str(i)])
    
    
%     klos_range = round(size(bcone,1)/2)+[-30:3:30];%round(size(bcone,1)/3):round(2*size(bcone,1)/3);
%     inspect_Erange = round((1/e_bin)*(250:500));
% %     figure, 
% %     subplot(211), pcolor(bcone_mod), shading flat
% %     subplot(212), plot(sum(bcone_mod,1))
%     ksym_search = zeros(1,length(klos_range));
%     for kli = 1:length(klos_range) 
%         klos = klos_range(kli);
%         kaxis_pre_sym = klos + round((1/k_bin)*(-100:100));    
%         [sym_bcone, sym_kaxis] = Symmetrized_spectra(bcone_mod(kaxis_pre_sym,:),round((1/k_bin)*(-100:100))');%1:size(cone_mod,1)]');
%         trial_dos = sum(repmat(abs(sym_kaxis),[1,length(inspect_Erange)]).*sym_bcone(:,inspect_Erange),1);
%         
% %         figure, 
% %         subplot(311), pcolor(bcone_mod), shading flat
% %         subplot(312), pcolor(sym_cone), shading flat
% %         subplot(313), plot(trial_dos)
%         ksym_search(kli) = min(trial_dos);
%     end
%     figure, plot(ksym_search)
%     the_klos = klos_range(ksym_search==min(ksym_search));
%     kaxis_pre_sym = the_klos + round((1/k_bin)*(-100:100));    

%     [sym_bcone,sym_kaxis] = Symmetrized_spectra(bcone_mod(kaxis_pre_sym,:),round((1/k_bin)*(-100:100))');
%     dos = sum(repmat(abs(sym_kaxis),[1,size(sym_bcone,2)]).*sym_bcone,1);

    %Finding average k
%     
%     avg_k = sum(sum(repmat([1:size(cone,1)]',[1,size(cone,2)]).*cone)) / sum(sum(cone));
%     
%     %Plotting DOS
%     binned_cone=Binning_2d(cone,dos_bin,1);
%     
%     intensity = sum(sum(binned_cone));
%     
%     int_k_range=50;
%     
%     dos=[(dos_bin*(0:floor(size(cone,2)/dos_bin)-1)+round(dos_bin/2))',(abs(-int_k_range:int_k_range)*binned_cone(round(avg_k)-int_k_range:round(avg_k)+int_k_range,:))'];
%     
%     %plot(dos(:,1),dos(:,2));
%     %pause(0.2);
%     
%     %%%------ fitting  start here (needs improvement) ---------%%%%%
%     
%     v_shape=1;
%     if v_shape==1   %fitting vshape dirac cone
%         int_E_range=[0,400];
%         dos_int=dos(dos(:,1)>int_E_range(1) & dos(:,1)<int_E_range(2),:);
%         
%         lower_cone=dos_int(:,1) < 150 & dos_int(:,1) > 0;
%         upper_cone=dos_int(:,1) > 200 & dos_int(:,1) < 400;
%         
%         %fit upper and lower cone
%         p1=fit(dos_int(lower_cone,1),dos_int(lower_cone,2),'poly1');
%         p2=fit(dos_int(upper_cone,1),dos_int(upper_cone,2),'poly1');
%         
%         if plot_fig==1
%             plot(p1,dos_int(:,1),dos_int(:,2))
%             hold on,
%             plot(p2,dos_int(:,1),dos_int(:,2));
%             hold off
%             legend('off');
%             pause(0.3);
%         end
%         
%         p1_val=coeffvalues(p1);
%         p2_val=coeffvalues(p2);
%         % calculate intersection
%         Dirac_energy=(p2_val(2)-p1_val(2))/(p1_val(1)-p2_val(1));
%         
%     else           %fitting straight and linear
%         int_E_range=[0,400];
%         dos_int=dos(dos(:,1)>int_E_range(1) & dos(:,1)<int_E_range(2),:);
%         
%         b_start=1;
%         dirac_start=175;
%         a_val=mean(dos_int(dos_int(:,1) < 150 & dos_int(:,1)>50,2));
%         
%         ft=fittype(['Dirac_line(x,',num2str(a_val),',b,E_dirac)']);
%         f=fit(dos_int(:,1),dos_int(:,2),ft,'Start',[b_start,dirac_start]);
%         f_coeff=coeffvalues(f);
%         Dirac_energy=f_coeff(1);
%         
%         if plot_fig==1
%             %figure();
%             plot(f,dos(:,1),dos(:,2));
%             text(250,6000,['DP = ', num2str(Dirac_energy)]);
%             pause(0.2);
%         end
%        
%     end
%     
%     %determine the momentum center near the Dirac point
%     if Dirac_energy>50 && Dirac_energy<size(cone,2)-10
%         momentum_curve=sum(binned_cone(:,floor(Dirac_energy/dos_bin)-2:floor(Dirac_energy/dos_bin)+2),2);
%         
%         k_range=(round(avg_k)-int_k_range:round(avg_k)+int_k_range)';
%         k_center=sum(k_range.*momentum_curve(k_range))/sum(momentum_curve(k_range));
%         
%         dirac_map(i)=Dirac_energy;
%         kc_map(i)=k_center;
%     end
    
    %%%%%-------------- fitting ends here ---------------- %%%%%%%
%end
%

%{
dirac_map = reshape(dirac_map,[X_size,Y_size]);

kc_map = reshape(kc_map,[X_size,Y_size]);

X=reshape(X,[X_size,Y_size]);
Y=reshape(Y,[X_size,Y_size]);

figure, imagesc(Y(1,:),X(:,1),dirac_map);
xlabel('Y (mm)');
ylabel('X (mm)');

%% adding each arpes image to get a full arpes image
%interpolated arpes images range
[E_interp,K_interp]=meshgrid(-200:0.5:250,-80:0.5:80);

%[E_interp,K_interp]=meshgrid(50:0.5:450,50:0.5:250);
ll=1;

full_arpes=zeros(size(E_interp));
anti_spectra=zeros(size(E_interp));

full_cone=zeros(size(cone));

sorted_dirac=sort(dirac_map(:));
edge = [sorted_dirac(1:floor(length(sorted_dirac)/3):end);sorted_dirac(end)];
n=0;

for jj= 1:Y_size
    for ii=1:X_size
        if kc_map(ii,jj)>0 && dirac_map(ii,jj)>0
            frame = reshape(data{1,end}(ii+(jj-1)*X_size,:), [768,997]);
            cone = frame(cone_range_K(1):cone_range_K(2), cone_range_E(1):cone_range_E(2));
            
            [E_coor,K_coor]=meshgrid(1:size(cone,2),1:size(cone,1));
            
            E_coor=E_coor-dirac_map(ii,jj);
            K_coor=K_coor-kc_map(ii,jj);
            
            interp_arpes=interp2(E_coor,K_coor,cone,E_interp,K_interp);
            %determine wheter to add this point into the sum or not
            if (abs(dirac_map(ii,jj)-mean(dirac_map(Neighbor_sites(ii,jj,X_size,Y_size))))<std(dirac_map(:)))...
                    %&& (dirac_map(ii,jj)>240) && dirac_map(ii,jj) <320
                full_arpes=full_arpes+interp_arpes;
            else
               n=n+1;
            end
            full_cone=full_cone+cone;
        end
    end
end

E_bin=10;
K_bin=3;
full_arpes_binned=Binning_2d(full_arpes,E_bin,K_bin);
E_binned = E_interp(1,1:E_bin:end-1);
if K_bin==1
    K_binned=K_interp(:,1);
else    
    K_binned = K_interp(1:K_bin:end-1,1);
end
full_dos_binned=abs(K_binned)'*full_arpes_binned;

figure,imagesc(E_binned,K_binned,full_arpes_binned.^0.3);axis xy;
%caxis([1000,5000])

figure, plot(E_binned,full_dos_binned);
%hold on
%}

%{
%% symmetrized binned ARPES spectra
sum_list=[];
spectra_list={};
binned_spectra=reshape(full_arpes_binned(full_arpes_binned>0),size(full_arpes_binned,1),[]);
for kpl=round(size(binned_spectra,1)/2)-5:round(size(binned_spectra,1)/2)+5
    asymmetrized_spectra=zeros(min([abs(kpl),abs(size(binned_spectra,1)-kpl)]),size(binned_spectra,1));
   
    asymmetrized_spectra=binned_spectra(kpl+1:kpl+size(asymmetrized_spectra,1),:)-flipud(binned_spectra(kpl-size(asymmetrized_spectra,1)+1:kpl,:));
    
    sym_spectra=(binned_spectra(kpl+1:kpl+size(asymmetrized_spectra,1),:)+flipud(binned_spectra(kpl-size(asymmetrized_spectra,1)+1:kpl,:)))/2;
    
    pcolor(asymmetrized_spectra),shading flat;
    pause(0.5);
    
    sum_list(end+1,:)=[kpl,sum(sum(asymmetrized_spectra(:).^2))];
    spectra_list(end+1)={sym_spectra};
end    

[m_val,m_kpl]=min(sum_list(:,2));
k_zero=sum_list(m_kpl,1);
full_sym_spectra=[flipud(spectra_list{m_kpl}(1:end-1,:));spectra_list{m_kpl}];
k_length=(size(full_sym_spectra,1)-1)/2;

Eaxis=E_binned(full_arpes_binned(1,:)>0);
Kaxis=[K_binned(k_zero-k_length:k_zero-1)-K_binned(k_zero);-flipud(K_binned(k_zero-k_length:k_zero)-K_binned(k_zero))];

sym_dos=abs(Kaxis)'*full_sym_spectra;

figure,imagesc(Eaxis,Kaxis,full_sym_spectra);axis xy;

figure, plot(Eaxis,sym_dos);

%%
%{
E_normalized=205;
figure, plot(E_binned,full_dos_binned/full_dos_binned(E_binned==200));
hold on;

plot(S_dos(:,1),S_dos(:,2)/S_dos(S_dos(:,1)==E_normalized,2))

plot(P_dos(:,1),P_dos(:,2)/P_dos(P_dos(:,1)==E_normalized,2)),
hold off;
toc;
%}

%% standard deviation of dirac enegy map
%{
zero_list=find(dirac_map==0);
dirac_nonzero=dirac_map;

for zpl=1:length(zero_list)
    n_sites=[zero_list(zpl)-X_size-1:zero_list(zpl)-X_size+1,zero_list(zpl)-1,zero_list(zpl)+1 ...
        zero_list(zpl)+X_size-1:zero_list(zpl)+X_size+1];
    
    dirac_nonzero(zero_list(zpl))=mean(dirac_map(dirac_map(n_sites)~=0));
end

dirac_map=dirac_nonzero;

std_list=[];
for rpl=1:size(dirac_map,1)
    std_list(end+1)=1/sqrt(2)*std(dirac_map(rpl,2:end)-dirac_map(rpl,1:end-1));
end


plot(1:31,std_list);
xlabel('Row number');
ylabel('Standard deviation');
%}
% measure deviations
%{
        deviation=mean(abs(dos_int(:,2)-feval(f,dos_int(:,1))));
        bad_list=find(abs(dos_int(:,2)-feval(f,dos_int(:,1)))>2*deviation);
        
        aa_val=mean(dos_int(dos_int(:,1) < 150 & ~ismember(dos_int(:,1),bad_list),2));
        
        ft_new=fittype(['Dirac_line(x,',num2str(aa_val),',b,E_dirac)']);
        f_new=fit(dos_int(:,1),dos_int(:,2),ft_new,'Start',[b_start,dirac_start],'Exclude',bad_list);
        f_coeff_new=coeffvalues(f_new);
         
        Dirac_energy=f_coeff_new(1);
        
        if plot_fig==1
            plot(f_new,dos(:,1),dos(:,2),'-');
            hold on, plot(dos_int(bad_list,1),dos_int(bad_list,2),'r*');
            text(100,max(dos_int(:,2))*2/3,['DP = ', num2str(Dirac_energy),char(10),'DP new = ',num2str(Dirac_new)]);
            hold off;
            pause(0.4);
        end
        %}
%}