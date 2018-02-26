if load_data == 1  
    load rfc_big_scan_170927.mat;
    load rfc_FL_scan_170927.mat;
    load rfc_ncorr_scan_171019w.mat;
    load cones.mat;
    load pre_dos_Es;
    load kLOS_171122.mat;
end

table_title = 'pre\_dos\_E fit';
dirac_Es = pre_dos_Es;%rfc_Es_after;   %input Energies vector (cone pixels)
dirac_ks = kLOS;   %input mtm vector (cone pixels)
EV_map = reshape(((rfc_FL_Es - dirac_Es) * pix2eV),31,31);


DPEs_eV = (rfc_FL_Es - pre_dos_Es)*pix2eV;%the_Es;%rfc_Es_after;%(rfc_FL_Es - rfc_small_cut_Es)*pix2eV;
B_map = reshape(DPEs_eV,31,31);   
%B_map = reshape(rfc_small_Es,31,31);

%E_interval_list = [350,360; 360,365; 365,370; 370,375; 375,380; 380,390; 390,410];%[370,385; 385,390; 390,395; 395,400; 400,405; 405,415];%; 386.7,396.9; 396.9,407.1; 407.1,422.4];% ; 351.5,356.6; 356.6,361.7; 361.7,366.8; 366.8,382];%,350; 350,355; 355,357.5; 357.5,360; 360,362.5; 362.5,365; 365,375; 375,380; 380,385 ];%[ .4,.43; .43,.45; .45,.46; .46,.47; .47,.49; .49,52 ];
E_interval_list = [0.425,0.45;0.45,0.46; 0.46,0.465; 0.465,0.47; 0.47,0.48; 0.48,0.49; 0.49,0.51];%[320,340; 340,345; 345,350; 350,352.5; 352.5,355; 355,360; 360,380];
B_interval_list = (E_interval_list - min(B_map(:)))/(max(B_map(:))-min(B_map(:)));

pre_filter = cat(1,[]);
pre_filter(1,:) = [174,189, kLOS];
pre_filter(2,:) = [220,440, DPI_big];
%pre_filter(3,:) = [9,27, rfc_corrspreads_after];
%pre_filter(3,:) = [.8,1, rfc_corrs_after];

figure
allfiltered = ones(1,961);% reshape(B_map,1,961);
for NN = 1:size(pre_filter,1)
    allfiltered(pre_filter(NN,3:end) < pre_filter(NN,1)) = 0;
    allfiltered(pre_filter(NN,3:end) > pre_filter(NN,2)) = 0;
    
    subplot(size(pre_filter,1),2,2*(NN-1)+2)
    histogram(pre_filter(NN,3:end),40), hold on;
    plot([pre_filter(NN,1),pre_filter(NN,1)],[0,max(hist(pre_filter(NN,3:end),40))],'r'), hold on;
    plot([pre_filter(NN,2),pre_filter(NN,2)],[0,max(hist(pre_filter(NN,3:end),40))],'r'), hold off;
end
filtered_A_map = reshape(allfiltered,31,31);
subplot(NN,2,[1]) 
imagesc(filtered_A_map), axis xy
A_map = filtered_A_map;
A_interval_list = [0.5,1];
%{
A_map = rfc_small_cut_ks;
A_interval_list = [173,176; 176-181; 182,187];
A_interval_list = (A_interval_list - min(A_map(:)))/(max(A_map(:))-min(A_map(:)));
%}

pix2eV = (1.599/(2*496));
pix2invA = 0.512*0.04631/180*3.1415*14/30*sqrt(110-4);

[E_interp,K_interp] = meshgrid(-199.5:.5:350, -80:.5:80); 
full_arpes = zeros(size(E_interp));

[region_list]=Correlation_list_indep_ranges(A_map,B_map,A_interval_list,B_interval_list);
regional_arpes=zeros([size(E_interp),size(region_list,1)]);

mean_DPE_eV = zeros(1,size(region_list,1));
mean_FL_eV = zeros(1,size(region_list,1));
map_regions = cat(3,[]);
NNN = 0;
for iii=1:size(region_list,1)
    nnn = 0;
    
     for jjj=find(region_list{iii}==1)'
         
         nnn = nnn+1;
         [map_y, map_x] = ind2sub([31,31],jjj);
         map_regions(1,nnn,iii) = map_x;
         map_regions(2,nnn,iii) = map_y;
         cone = result1i(:,:,jjj);
         [E_coor, K_coor] = meshgrid(1:size(cone,2),1:size(cone,1));
         E_coor = E_coor - rfc_FL_Es(jjj);%dirac_Es(jjj);  
         K_coor = K_coor - dirac_ks(jjj);  
         interp_arpes = interp2(E_coor, K_coor, cone, E_interp, K_interp);
         regional_arpes(:,:,iii) = regional_arpes(:,:,iii) + interp_arpes;
         
         DPE_pix = dirac_Es(jjj);
         FL_pix = rfc_FL_Es(jjj);
         DPE_eV = (-FL_pix + DPE_pix) * pix2eV;
         FL_eV = FL_pix * pix2eV;
         mean_DPE_eV(iii) = mean_DPE_eV(iii) + DPE_eV;
         mean_FL_eV(iii) = mean_FL_eV(iii) + FL_eV;
     end
  
         
     mean_DPE_eV(iii) = mean_DPE_eV(iii) / length(find(region_list{iii}==1));
     mean_FL_eV(iii) = mean_FL_eV(iii) / length(find(region_list{iii}==1));
     nnn_scans(iii) = nnn;
     NNN = NNN + nnn;
end


% plot out the binned and symmetrized regional spectra

E_bin=10;
K_bin=3;

E_binned = E_interp(1,1:E_bin:end-1);

if K_bin == 1
    K_binned = K_interp(:,1);
else    
    K_binned = K_interp(1:K_bin:end,1);
end
Eaxis=E_binned;
arp_chart = figure;
nn=0;
nnn=0;
for iii=1:size(region_list,1)    
        
    region_arpes_binned = Binning_2d(regional_arpes(:,:,iii),E_bin,K_bin);
    [region_arpes_symmetrized,Kaxis]=Symmetrized_spectra(region_arpes_binned,K_binned);
    ax = subplot(2,size(B_interval_list,1),iii); %%%%%%replace 2 with size(B_interval_list,1)%%%
    
    %imagesc(regional_arpes(:,:,iii)), axis xy
    %imagesc(rot90(Eaxis,-1),rot90(Kaxis,-1),rot90(region_arpes_binned(1:103,:),-1)), %axis xy;
    %imagesc(Eaxis,Kaxis,region_arpes_symmetrized), axis xy;
    Eaxis_eV = E_interp(1,1:E_bin:end-1) * pix2eV;
    Eaxis_eV_0fixed = Eaxis_eV - abs(mean_DPE_eV(iii));           
    Kaxis_invA = Kaxis * pix2invA;
    if length(find(region_list{iii}==1)) > 0
        %imagesc(Kaxis_invA, flip(Eaxis_eV_0fixed), rot90(norman(region_arpes_symmetrized,0,5))), axis xy
        imagesc(Kaxis_invA, flip(Eaxis_eV_0fixed), rot90(norman(region_arpes_binned(1:103,:),0,5))), axis xy

        
        %yticks([0,-(mean_DPE_eV(iii)/pix2eV)])
        %yticklabels({num2str(mean_DPE_eV(iii)),'0'})
        ylabel_str = 'E - E$_{F}$, (eV)';
        ylabel(ylabel_str,'Interpreter','latex','FontSize',8)
        %yt = get(gca, 'YTick');
        %set(gca, 'FontSize',7)
        %ytickangle(90)
        ylim([-0.75,0.1])
        %xticks([-50,0,50])
        %xticklabels({num2str(-50*pix2invA),'0',num2str(50*pix2invA)})
        xlabel('Momentum ({\AA}$^{-1}$)','Interpreter','latex','FontSize',8)
        %xt = get(gca, 'XTick');
        %set(gca, 'FontSize', 7)
    end
    if iii == 1+size(B_interval_list,1)*nn
        descrr = {[num2str(A_interval_list(nn+1,1)),'-',num2str(A_interval_list(nn+1,2))]};
        axes(ax)
        text(-0.2,-0.3,descrr')
        %ylabel([num2str(A_interval_list(nn+1,1)),'-',num2str(A_interval_list(nn+1,2))]);
        nn=nn+1;
    end 
    if ismember(iii,[1:size(B_interval_list,1)]) == 1
        %set(gca,'XaxisLocation','top');
        descr = {[num2str(E_interval_list(nnn+1,1)),'-',num2str(E_interval_list(nnn+1,2))]};
        axes(ax)
        text(-.05,.15,descr)
        %xlabel([num2str(B_interval_list(nnn+1,1)),'-',num2str(B_interval_list(nnn+1,2))]);
        nnn=nnn+1;
    end
    ax2 = subplot(2,size(B_interval_list,1),iii+size(B_interval_list,1));
    imagesc(EV_map), axis xy, hold on;
    plot(map_regions(1,:,iii),map_regions(2,:,iii),'k*'), hold off;
    title(['nnn=',num2str(nnn_scans(iii))])
    colormap jet
    caxis(ax2, [.38,.65])
end
suptitle([table_title,', total scans: ',num2str(NNN)]);




