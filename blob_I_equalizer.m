
function [blob_I_normed_cones, pix_positions,ave_blob_intensities] = blob_I_equalizer(frame_data, frame_size, num_scans, cone_range_K, cone_range_E, box_row_range, box_col_width, box_col_lowlim, box_col_uplim) 
tic;
blob_I_normed_cones = zeros(length(cone_range_K(1):cone_range_K(end)),length(cone_range_E(1):cone_range_E(end)),num_scans);
not_normed_cones = zeros(size(blob_I_normed_cones));

box_col_ranges = {};
box_col_start = box_col_lowlim;
box_col_end = box_col_start + box_col_width - 1;
box_col_i = 1;
while box_col_end <= box_col_uplim
    box_col_ranges{box_col_i} = (box_col_start:1:box_col_end);
    box_col_i = box_col_i + 1;
    box_col_start = box_col_start + box_col_width;
    box_col_end = box_col_start + box_col_width - 1;
end
box_col_ranges = box_col_ranges';

ave_blob_intensities = zeros(1,size(box_col_ranges,1));
pix_positions = zeros(1,size(box_col_ranges,1));

for fcr_i = 1:size(box_col_ranges,1)
    if rem(fcr_i,ceil(size(box_col_ranges,1)/4)) == 0
        disp(['Analyzing frame blobs ',num2str( round((fcr_i/(size(box_col_ranges,1)))*100)),'% done'])
    end    
    box_col_range = box_col_ranges{fcr_i};
    
    numbers_of_blobs = NaN*ones(1,num_scans);

    all_blobs = cell(1,num_scans);
    all_blob_sizes = zeros(1,100000000);
    all_blob_Is = zeros(1,100000000);
    all_blob_i = 1;
    
    for i = 1:num_scans
        frame = reshape(frame_data{1,end}(i,:), frame_size);
        not_normed_cones(:,:,i) = frame(cone_range_K(1):cone_range_K(end),cone_range_E(1):cone_range_E(end));
        
        box = frame(box_row_range,box_col_range);
        box01 = box;
        box01(box01>0) = 1;
        blob = bwconncomp(box01,4);
        number_blobs = blob.NumObjects;
        numbers_of_blobs(i) = number_blobs;
        blob_specs = zeros(2,number_blobs); %1. no pixels, 2.total intensity

        for blob_i = 1:number_blobs
            blob_in_question = blob.PixelIdxList{blob_i};
            number_pixels = length(blob_in_question);
            pixel_intensity = zeros(number_pixels,1);
            blob_intensity = 0;
            for blob_ii = 1:number_pixels
                blob_pix_ind = blob_in_question(blob_ii);
                [blob_pix_row, blob_pix_col] = ind2sub(size(box), blob_pix_ind);
                pixel_intensity(blob_ii) = box(blob_pix_row,blob_pix_col);        
                blob_intensity = blob_intensity + pixel_intensity(blob_ii);
            end     
            blob_specs(1,blob_i) = number_pixels;
            blob_specs(2,blob_i) = blob_intensity;

        end
        all_blobs{i} = blob_specs; 
        all_blob_Is(all_blob_i:all_blob_i+number_blobs-1) = blob_specs(2,:);
        all_blob_sizes(all_blob_i:all_blob_i+number_blobs-1) = blob_specs(1,:);
        all_blob_i = all_blob_i + number_blobs;
    end
    ave_blob_intensity = mean(all_blob_Is(all_blob_Is>0));
    
    ave_blob_intensities(fcr_i) = ave_blob_intensity;
    pix_positions(fcr_i) = round(median(box_col_range));
    
end

interped_ave_blob_Is = interp1(pix_positions,ave_blob_intensities,[cone_range_E(1):cone_range_E(end)],'pchip');
blob_I_normalizer_mat = repmat(interped_ave_blob_Is,[length(cone_range_K(1):cone_range_K(end)),1]);

disp('Blob intensities analyzed, starting cone blob I normalization')
for i = 1:num_scans   
    blob_I_normed_cones(:,:,i) = not_normed_cones(:,:,i)./blob_I_normalizer_mat;
end
toc;
end

