function [ neighbour_list ] = Neighbor_sites( ii,jj,X_size,Y_size )
%NEIGHBOUR_SITES Summary of this function goes here
%   Detailed explanation goes here
    
    [n_x,n_y]=meshgrid(ii-1:ii+1,jj-1:jj+1);
    neighbour_list=[reshape(n_x,[numel(n_x),1]),reshape(n_y,[numel(n_y),1])];
    neighbour_list=neighbour_list(~ismember(neighbour_list,[ii,jj],'rows'),:);
    
    good_list=(neighbour_list(:,1)>=1 & neighbour_list(:,1) <=X_size) .* (neighbour_list(:,2)>=1 & neighbour_list(:,2) <=Y_size);
    neighbour_list=neighbour_list(logical(good_list),:);
    
    neighbour_list=neighbour_list(:,1)+(neighbour_list(:,2)-1)*X_size;
    
end

