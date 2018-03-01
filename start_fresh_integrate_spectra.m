load_data = 1;
if load_data == 1  
    load 'start_fresh_cluster1_180228.mat';
    load 'normed_binned_cones.mat';
    load 'rfc_FL_scan_170927.mat';
    load 'top_spans.mat';
    load 'mins_bel_bcb.mat';
    load 'cool_purply_science_colormap.mat';
end

num_scans = 961; e_bin = 3; k_bin = 1;
AVE_FL = nanmean(rfc_FL_Es);
pix2eV = (1.599/(2*496));
pix2invA = 0.512*0.04631/180*3.1415*14/30*sqrt(110-4);
table_title = 'Cluster1 Scan 180228';
dirac_Es = dirac_E_map;%bcb_finds;%BCB_Es;%rfc_Es_after;   %input Energies vector (cone pixels)
dirac_ks = dirac_K_map;%kLOS;%.5*(rfc_ks_after+kLOS);   %input mtm vector (cone pixels)

DPEs_eV = (mean(rfc_FL_Es) - dirac_Es)*pix2eV;%(rfc_FL_Es - bcb_findss)*pix2eV;%bcb_finds;%rfc_Es_after;%(rfc_FL_Es - rfc_small_cut_Es)*pix2eV;
B_map = reshape(DPEs_eV,31,31);   
%B_map = reshape(dirac_V_map*pix2eV/pix2invA,31,31);

B_interval_list_abs = [.42,.44; .44,.46; .46,.48; .48,.50; .50,.52];
B_interval_list_abs = [.40,.42; .42,.44; .44,.46; .46,.48; .48,.50; .50,.52; .52,.54];%; .445,.455; .455,.465; .465,.475; .475, .495; .495,.505]; 
%B_interval_list_abs = [2.8,3; 3.0,3.2; 3.2,3.4; 3.4,3.6; 3.6,3.8; 3.8,4.0];%3.1,3.4; 3.4,3.7; 3.7,4.0; 4.0,4.3];
B_interval_list = (B_interval_list_abs - min(B_map(:)))/(max(B_map(:))-min(B_map(:)));

pre_filter = cat(1,[]);
pre_filter(1,:) = [171,189, dirac_ks];
pre_filter(2,:) = [200,inf, DP_intensity_map];
%pre_filter(3,:) = [.15,inf, fit_evals_map];
%pre_filter(3,:) = [0,5, eval_under_bcb];
%pre_filter(4,:) = [mean(top_spans)-std(top_spans), mean(top_spans)+std(top_spans), top_spans];
pre_filter(3,:) = [0,1 mins_bel_bcb];
pre_filter(4,:) = [.2,inf, fit_evals_map];

figure
allfiltered = ones(1,961);%reshape(B__map,1,961);
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

[binE_interp,binK_interp] = meshgrid(-round(600/e_bin):1:round(80/e_bin), -round(100/k_bin):1:round(100/k_bin)); 
full_arpes = zeros(size(binE_interp));

[region_list]=Correlation_list_indep_ranges(A_map,B_map,A_interval_list,B_interval_list);
number_panels = size(region_list,1);

regional_arpes=zeros([size(binE_interp),number_panels]);

interped_arpes = zeros([size(binE_interp),num_scans]);
involved_scans = zeros(1,961);
involved_scans_panel = zeros(1,961,number_panels);
mean_DPE_eV = zeros(1,number_panels);
mean_FL_eV = zeros(1,number_panels);
map_regions = cat(3,[]);
NNN = 0;
panel_nnn = zeros(1,number_panels);
scan_box_events = NaN*ones(1,961,number_panels);
panel_mean = zeros(1,number_panels);
for iii=1:number_panels
    nnn = 0;
    
     for jjj=find(region_list{iii}==1)'
         involved_scans(jjj) = involved_scans(jjj)+1;
         involved_scans_panel(1,jjj,iii) = involved_scans_panel(1,jjj,iii)+1;
         nnn = nnn+1;
         [map_y, map_x] = ind2sub([31,31],jjj);
         map_regions(1,nnn,iii) = map_x;
         map_regions(2,nnn,iii) = map_y;
         
         bcone = binned_cones(:,:,jjj);
         [bE_coor, bK_coor] = meshgrid(1:size(bcone,2),1:size(bcone,1));
         bE_coor = bE_coor - mean(rfc_FL_Es)/e_bin;%rfc_FL_Es(jjj)/e_bin;%dirac_Es(jjj)/e_bin;  
         bK_coor = bK_coor - dirac_ks(jjj)/k_bin;  
         interp_arpes = interp2(bE_coor, bK_coor, bcone, binE_interp, binK_interp);
         interped_arpes(:,:,jjj) = interp_arpes;
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
     panel_nnn(iii) = nnn;
     NNN = NNN + nnn;
end
%%

arp_chart = figure;
nn=0;
nnn=0;
out_spec = cell([1,size(region_list,1)]);
out_kax = cell([1,size(region_list,1)]);
out_eax = cell([1,size(region_list,1)]);
out_nnn = zeros(1,size(region_list,1));

for iii=1:size(region_list,1)    
    regional_arpes(:,:,iii) = imgaussfilt(regional_arpes(:,:,iii),1); 
    region_arpes = regional_arpes(:,:,iii);
    [region_arpes_symmetrized, sym_Kaxis]=Symmetrized_spectra(region_arpes,binK_interp(:,1));    
    
    Eaxis_eV_0fixed = binE_interp(1,:) * e_bin * pix2eV;
    Kaxis_invA = binK_interp(:,1) * k_bin * pix2invA;    
    sym_Kaxis_invA = sym_Kaxis * k_bin * pix2invA;
    
    ax = subplot(2,size(B_interval_list,1),iii); 
    if length(find(region_list{iii}==1)) > 0
        
        %imagesc(sym_Kaxis_invA, flip(Eaxis_eV_0fixed), rot90(norman(region_arpes_symmetrized,0,10))), axis xy
        imagesc(Kaxis_invA, flip(Eaxis_eV_0fixed), rot90(norman(region_arpes,0,30))), axis xy
        
        %yticks([0,-(mean_DPE_eV(iii)/pix2eV)])
        %yticklabels({num2str(mean_DPE_eV(iii)),'0'})
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
        %descrr = {[num2str(A_interval_list(nn+1,1)),'-',num2str(A_interval_list(nn+1,2))]};
        axes(ax)
        ylabel_str = 'E - E$_{F}$, (eV)';
        ylabel(ylabel_str,'Interpreter','latex','FontSize',8)
        %text(-0.2,-0.3,descrr')
        %ylabel([num2str(A_interval_list(nn+1,1)),'-',num2str(A_interval_list(nn+1,2))]);
        nn=nn+1;
    end 
    
    if ismember(iii,[1:size(B_interval_list,1)]) == 1
        %set(gca,'XaxisLocation','top');
        descr = {[num2str(B_interval_list_abs(nnn+1,1)),'-',num2str(B_interval_list_abs(nnn+1,2))]};
        axes(ax)
        text(-.05,.15,descr)
        %xlabel([num2str(B_interval_list(nnn+1,1)),'-',num2str(B_interval_list(nnn+1,2))]);
        nnn=nnn+1;
    end
    
    
    map_of_color = [0,0,0;43,15,107;93,0,200;198,0,116;235,83,60;245,151,48;233,216,57;255,255,255]/255;
       
    ax2 = subplot(2,size(B_interval_list,1),iii+size(B_interval_list,1));
    imagesc(B_map), axis xy, hold on; %caxis([min(B_map(B_map>0)),max(B_map(B_map>0))]), hold on;
    plot(map_regions(1,:,iii),map_regions(2,:,iii),'p','Color',[0,255,255]/255,'LineWidth',1), hold off;
    title(['nnn=',num2str(panel_nnn(iii))])
    colormap(map_of_color_64)
    
    output_kaxis_invA = Kaxis_invA;
    
    out_spec{1,iii} = region_arpes;
    out_kax{1,iii} = output_kaxis_invA;
    out_eax{1,iii} = Eaxis_eV_0fixed;
    out_nnn(iii) = panel_nnn(iii);
end
suptitle([table_title,', total scans: ',num2str(NNN)]);





