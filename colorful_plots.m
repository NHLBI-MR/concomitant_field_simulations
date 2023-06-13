function [h,lgnd] = colorful_plots(x,y,lgnd)
% [h,lgnd] = colorful_plots(x,y <opt.>,lgnd <opt.>)
% Line plots arranged in a parula order (padded to avoid yellow on white
% lines)
% Auto legend can be removed 

show_legend = 0;
if nargin == 1
    y = x;
    x = repmat(1:size(y,1),[size(y,2),1])';
elseif ((nargin==2) && iscell(y))
    lgnd = y;
    show_legend = 1;
    y = x;
    x = repmat(1:size(y,1),[size(y,2),1])';
elseif (nargin==3)
    show_legend = 1;
    if ~iscell(lgnd)
        % vector input
        lgnd = cellfun(@num2str,num2cell(lgnd), 'UniformOutput', 0);
    end
end

% h = plot(x, y, '-');
h = plot(x, y, '-', 'LineWidth', 1.5);
% axis image;

% pad parula to avoid yellows
p_cmap = parula(length(h)+3); clear clz;

for i=1:length(h)
    clz{i,1} = p_cmap(i,:);
    
    % auto legend 
    if nargin < 3 && ~show_legend
       lgnd{i} = num2str(i); 
    end
end

set(h, {'color'}, clz);

if show_legend
legend(lgnd, 'Location', 'best');
end

end