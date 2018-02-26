function [ region_list ] = Correlation_list_indep_ranges( A_map,B_map,A_ranges,B_ranges)
%CORRELATION_PLOT Summary of this function goes here
%   Detailed explanation goes here split points should be a n by 2 matrix
%   indentifying the start and end point for each interval 
%   A_map and B_map are the two maps, output region_list is a n*n by 1 cell
%   containing the points of each region
region_list=cell(size(A_ranges,1)*size(B_ranges,1),1);
A_map_norm=(A_map-min(A_map(:)))/(max(A_map(:))-min(A_map(:)));
B_map_norm=(B_map-min(B_map(:)))/(max(B_map(:))-min(B_map(:)));
for ii=1:size(A_ranges,1)  % ii is for A_map
    for jj = 1:size(B_ranges,1) % jj is for B_map
        list_pl=(ii-1)*size(B_ranges,1)+jj;
        region_list{list_pl}= (A_map_norm >= A_ranges(ii,1) & A_map_norm <= A_ranges(ii,2)) .* (B_map_norm >= B_ranges(jj,1) & B_map_norm <= B_ranges(jj,2));
    end
end



end
