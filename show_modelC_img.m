function [ model ] = show_modelC_img( cls, varargin )

set_conf_img;
set_param_GM;

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'conf'
      conf = arg;  
    case 'verbose'
      verbose = arg;
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

load(getModelPath(conf,cls), 'model');

svm = model.svm;
model_o = init_svm_img(model);
svm_o = model_o.svm;

% speficy graph type
frame_m = model.view_m.frame;
color_m = model.colorCode;
adjmat_m = model.graph_m.adjmat;

%svm = renormA(svm);
wV = sqrt(sum(svm.wV.^2,1)); 
wE = sqrt(sum(svm.wE.^2,1)); 
%wV = sum(abs(svm.wV),1); 
%wE = sum(abs(svm.wE),1); 

% show node shapes (ellipses)
% for i = 1:size(frame_m,2)
%     fprintf('%f - node  %d (feat %d)\n',wV(i),i,model.view_m.type(i));
% end
% for k = 1:size(edges_m,2)
%     fprintf('%f %f %f - edge (%d,%d)\n',svm.aE(1:3,k),edges_m(1,k),edges_m(2,k));
% end
% fprintf('\n%f - bias \n',svm.aBias);

% show the image
figure(1); 
imshow(rgb2gray(model.view_m.img)); hold on;
drawWeightedGraph( frame_m, [ 0 0 ], color_m, adjmat_m, wV, wE(:) );

% show GT bbox
%box_handle = rectangle('position', [box_coord(3), box_coord(1), box_coord(4)-box_coord(3), box_coord(2)-box_coord(1)]);
%set(box_handle, 'edgecolor','y', 'linewidth',5);    

%visualizeSIFT(model.view_m.desc, frame_m);
%visualizeFeatures(frame_m, 'style', 'frame', 'colorcode', 'g');
%visualizeSIFT(svm.aV, frame_m);

featBin_o = svm_o.wE;
rho_bins = svm_o.rho_bins;
theta_bins = svm_o.theta_bins;
featBin = svm.wE;

% normalizing
featBin_o = featBin_o ./ repmat(sum(abs(featBin_o),1), [ size(featBin_o,1) 1]); 
featBin = featBin ./ repmat(sum(abs(featBin),1), [ size(featBin,1) 1]); 

rb_o = featBin_o(1:numel(rho_bins),:);
tb_o = featBin_o(numel(rho_bins)+1:numel(rho_bins)+numel(theta_bins),:);  
rb = featBin(1:numel(rho_bins),:);
tb = featBin(numel(rho_bins)+1:numel(rho_bins)+numel(theta_bins),:);

for i=1:numel(rho_bins)
    rho_labels{i} = sprintf('< %.1f',rho_bins(i));
end
for i=1:numel(theta_bins)
    theta_labels{i} = sprintf('%.0f',theta_bins(i)*180/pi);
end
%theta_labels = 1:numel(theta_bins);

[ valE indE] = sort(wE,'descend');
for k=1:numel(indE)
    [ i j ] = ind2sub(size(adjmat_m),indE(k));
    fprintf('wE: %f\n',valE(k));
%for i=1:size(frame_m,2)
%     d_ori = sqrt(abs(frame_m(3,i)*frame_m(6,i)-frame_m(4,i)*frame_m(5,i)));
%     for k = 1:numel(rho_bins)
%         drawCircle(frame_m(1,i),frame_m(2,i), rho_bins(k)*d_ori, 'r-', 'linewidth', 2);
%     end
%     x = frame_m(1,i);
%     y = frame_m(2,i);
%     r = rho_bins(end-1)*d_ori;
%     pline_x = r * cos(theta_bins) + x;
%     pline_y = r * sin(theta_bins) + y;
%     for k=1:numel(theta_bins)
%         plot( [x pline_x(k)], [y pline_y(k)], 'r', 'linewidth', 2);
%     end
    
    %for j=1:size(frame_m,2)
        figure(1); hold on; 
        plot(frame_m(1,i), frame_m(2,i), 'o','Color', 'r', 'MarkerSize', 10 );
        plot(frame_m(1,[i j]), frame_m(2,[i j]), '-','Color', 'r', 'LineWidth', 5 );
        
        figure(2);
        %fprintf('%f %f %f\n', 180*ang/pi, dist, exp(area*log(2)));
        subplot(5,1,1);
        cla; hold on;
        barRho_o = bar(rb_o(:,i+(j-1)*size(frame_m,2)), 'FaceColor', 'b', 'EdgeColor', 'b'); 
        set(barRho_o,'BarWidth',0.8);
        barRho = bar(rb(:,i+(j-1)*size(frame_m,2)), 'FaceColor', 'r', 'EdgeColor', 'r');
        set(barRho,'BarWidth',0.4); 
        set(gca,'XTick', 1:numel(rho_bins) );
        set(gca,'XTickLabel', rho_labels );
        hold off;
        
        subplot(5,1,2);
        cla; hold on;
        barTheta_o = bar(tb_o(:,i+(j-1)*size(frame_m,2)), 'FaceColor', 'b', 'EdgeColor', 'b');
        set(barTheta_o,'BarWidth',0.5);
        barTheta = bar(tb(:,i+(j-1)*size(frame_m,2)), 'FaceColor', 'r', 'EdgeColor', 'r');
        set(barTheta,'BarWidth',0.25);
        set(gca,'XTick', 1:numel(theta_bins) );
        set(gca,'XTickLabel', theta_labels );
        hold off;
        
        subplot(5,1, [3 4 5]);
        cla; hold on;
        roseBar(theta_bins, tb_o(:,i+(j-1)*size(frame_m,2)), tb(:,i+(j-1)*size(frame_m,2)) );
        hold off;
        
        pause;
        
        figure(1);
        plot(frame_m(1,i), frame_m(2,i), 'o','Color', 'k', 'MarkerSize', 10 );
        plot(frame_m(1,[i j]), frame_m(2,[i j]), '-','Color', 'k', 'LineWidth', 5 );
        
    %end
    %pause;
end


% show the image
%figure('Name','Model Graph', 'NumberTitle', 'off'); 
%imshow(rgb2gray(model.view_m.img)+150); hold on;

% visualize nodes with their weights
%drawWeightedGraph( frame_m, [ 0 0 ], color_m, edges_m, wV, wE );

hold off


end




