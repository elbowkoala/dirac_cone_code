tic;
name = '20170518_00035.fits';
%info = fitsinfo(name);

%data = fitsread(name, 'BinaryTable');
load('P pol data.mat');
%%
plot_fig=0;

X_size = 31;
Y_size = 31;
num_scans = X_size *Y_size;

dirac_new_map = zeros(1,num_scans);
dirac_map = zeros(1,num_scans);
kc_map = zeros(1, num_scans);

X = zeros(1, num_scans);
Y = zeros(1, num_scans);

dos_bin=5;
cone_range_K=[400,700];
cone_range_E=[200,650];


for i = 1:1:num_scans
    
    X(i) = data{1,2}(i);
    Y(i) = data{1,3}(i);
    
    %plotting cone
    frame = reshape(data{1,end}(i,:), [768,997]);
    cone = frame(cone_range_K(1):cone_range_K(2), cone_range_E(1):cone_range_E(2));
    if plot_fig == 1
        %figure,pcolor(cone),shading flat;
        %pause(0.2);
    end
    %Finding average k
    
    avg_k = sum(sum(repmat([1:size(cone,1)]',[1,size(cone,2)]).*cone)) / sum(sum(cone));
    
    %Plotting DOS
    binned_cone=Binning_2d(cone,dos_bin,1);
    
    intensity = sum(sum(binned_cone));
    
    int_k_range=50;
    
    dos=[(dos_bin*(0:floor(size(cone,2)/dos_bin)-1)+round(dos_bin/2))',(abs(-int_k_range:int_k_range)*binned_cone(round(avg_k)-int_k_range:round(avg_k)+int_k_range,:))'];
    
    %plot(dos(:,1),dos(:,2));
    %pause(0.2);
    
    %%%------ fitting  start here (needs improvement) ---------%%%%%
    
    v_shape=1;
    if v_shape==1   %fitting vshape dirac cone
        int_E_range=[0,400];
        dos_int=dos(dos(:,1)>int_E_range(1) & dos(:,1)<int_E_range(2),:);
        
        lower_cone=dos_int(:,1) < 150 & dos_int(:,1) > 0;
        upper_cone=dos_int(:,1) > 200 & dos_int(:,1) < 400;
        
        %fit upper and lower cone
        p1=fit(dos_int(lower_cone,1),dos_int(lower_cone,2),'poly1');
        p2=fit(dos_int(upper_cone,1),dos_int(upper_cone,2),'poly1');
        
        if plot_fig==1
            plot(p1,dos_int(:,1),dos_int(:,2))
            hold on,
            plot(p2,dos_int(:,1),dos_int(:,2));
            hold off
            legend('off');
            pause(0.3);
        end
        
        p1_val=coeffvalues(p1);
        p2_val=coeffvalues(p2);
        % calculate intersection
        Dirac_energy=(p2_val(2)-p1_val(2))/(p1_val(1)-p2_val(1));
        
    else           %fitting straight and linear
        int_E_range=[0,400];
        dos_int=dos(dos(:,1)>int_E_range(1) & dos(:,1)<int_E_range(2),:);
        
        b_start=1;
        dirac_start=175;
        a_val=mean(dos_int(dos_int(:,1) < 150 & dos_int(:,1)>50,2));
        
        ft=fittype(['Dirac_line(x,',num2str(a_val),',b,E_dirac)']);
        f=fit(dos_int(:,1),dos_int(:,2),ft,'Start',[b_start,dirac_start]);
        f_coeff=coeffvalues(f);
        Dirac_energy=f_coeff(1);
        
        if plot_fig==1
            %figure();
            plot(f,dos(:,1),dos(:,2));
            text(250,6000,['DP = ', num2str(Dirac_energy)]);
            pause(0.2);
        end
       
    end
    
    %determine the momentum center near the Dirac point
    if Dirac_energy>50 && Dirac_energy<size(cone,2)-10
        momentum_curve=sum(binned_cone(:,floor(Dirac_energy/dos_bin)-2:floor(Dirac_energy/dos_bin)+2),2);
        
        k_range=(round(avg_k)-int_k_range:round(avg_k)+int_k_range)';
        k_center=sum(k_range.*momentum_curve(k_range))/sum(momentum_curve(k_range));
        
        dirac_map(i)=Dirac_energy;
        kc_map(i)=k_center;
    end
    
    %%%%%-------------- fitting ends here ---------------- %%%%%%%
end
%
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
