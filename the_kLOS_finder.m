
function kLOS_pt = the_kLOS_finder( bcone, kLOS_sigma, kLOS_erange, I_low_thresh, I_up_thresh)

%kLOS_erange = round((1/e_bin)*(200:550));
%wannasee = 1;
%kLOS_sigma = 5;

%tic
%kLOS = zeros(1,961);
%I_up_thresh = 4;
%I_low_thresh = 0;


    
    %scone = binned_cones(:,:,i);
scone = bcone;
sconee = imgaussfilt(scone, kLOS_sigma);
sconeee = sconee;
sconeee(sconee<(mean(sconee(:))-I_low_thresh*std(sconee(:))))=0;
sconeee(sconee>(mean(sconee(:))+I_up_thresh*std(sconee(:)))) = mean(scone(:))+I_up_thresh*std(sconee(:));
%figure, imagesc(sconeee), axis xy

kLOS_krange = round(size(scone,1)/3 : 2*size(scone,1)/3);%(100:250);
normcorr = zeros(1,length(kLOS_krange));
krow_n = 1;
for krow =  kLOS_krange
    hwidth = min(krow-1, (size(sconee,1)-krow-1)); %60;
    f_img1 = mat2gray(sconee(krow-hwidth:krow+hwidth, kLOS_erange));
    f_img2 = flipud(f_img1);

    f_mean = mean(f_img1(:));
    f_denom = sum(sum((f_img1-f_mean).^2));

    normcorr(krow_n) = sum(sum((f_img1-f_mean).*(f_img2-f_mean)))/f_denom;
    krow_n = krow_n+1;
end
%figure, imagesc(normcorr), axis xy

kLOS_pt = kLOS_krange(1) + find(normcorr==max(normcorr)) -1;

% if wannasee == 1
%   %  disp(['i=',num2str(i),' kLOS=',num2str(kLOS_pt)])
% 
%     figure, 
%     subplot(2,2,1), imagesc([1:size(scone,1)],fliplr([1:size(scone,2)]),rot90(imgaussfilt(scone,3))), axis xy, hold on;
%     plot([kLOS_pt,kLOS_pt],[1,size(scone,2)],'r'), hold on;
%     %plot([kLOS_b4(i),kLOS_b4(i)],[1,size(scone,2)],'w');
%     subplot(2,2,2), imagesc([1:size(sconee,1)],fliplr([1:size(sconee,2)]),rot90(sconeee,1)), axis xy, hold on;
%     %plot([0,size(sconee,2)],[kLOS_pt,kLOS_pt],'r'), hold on;
%     rectangle('Position',[kLOS_pt-hwidth,kLOS_erange(1),2*hwidth,kLOS_erange(end)-kLOS_erange(1)],'LineStyle',':','EdgeColor','r'), hold on;
%     plot([kLOS_pt,kLOS_pt],[0,size(sconee,2)],'r'), hold on;
%     %plot([0,size(sconee,2)],[kLOS_b4(i),kLOS_b4(i)],'w'), hold off
%     %plot([kLOS_b4(i),kLOS_b4(i)],[0,size(sconee,2)],'w'), hold off
%     subplot(2,2,[3,4]), plot(normcorr)
%     %title(['i=',num2str(i),' kLOSb4=',num2str(kLOS_b4(i)),' kLOS=',num2str(kLOS_pt)])
% end

end

