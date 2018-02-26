
function [X_finalized] = X_editor3(XVFg_matrix,a3,a4,draw_y_range)
% XVFg_matrix = XVFg_matrices{1,1};
% a3=380;
% a4=174;
% draw_y_range = 334:522;

X_mat = XVFg_matrix;
matsize1 = size(XVFg_matrix,1);
matsize2 = size(XVFg_matrix,2);


if a3-round(matsize2/2) > 0
    X_mat = horzcat(zeros(matsize1,a3-round(matsize2/2)),X_mat(:,1:end-(a3-round(matsize2/2))));
elseif a3-round(matsize2/2) < 0
    X_mat = horzcat(X_mat(:,round(matsize2/2)-a3+1:end),zeros(matsize1,round(matsize2/2)-a3));
end

if a4-round(matsize1/2) > 0   
    X_mat = vertcat(zeros(a4-round(matsize1/2),matsize2),X_mat(1:end-(a4-round(matsize1/2)),:));
elseif a4-round(matsize1/2) < 0
    X_mat = vertcat(X_mat(round(matsize1/2)-a4+1:end,:),zeros(round(matsize1/2)-a4,matsize2));
end

% draw_x_range1 = find(X_mat(:,draw_y_range(end))>0,1,'first');
% draw_x_range2 = find(X_mat(:,draw_y_range(end))>0,1,'last');

%X_matt = X_mat(draw_x_range1:draw_x_range2,draw_y_range);

%crop_x_range = draw_x_range1:draw_x_range2;
%crop_y_range = draw_y_range;

X_finalized = X_mat(:,draw_y_range);



end


