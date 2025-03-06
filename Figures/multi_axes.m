function [x_st, y_st, x_w, y_h] = multi_axes(fw, fh, nc, nr, ratio, margins, wmargins)
% [x_st, y_st, x_w, y_h] = multi_axes(fw, fh, nc, nr, ratio, margins)
% margins: [l_margin, r_margin, t_margin, b_margin]
% get axes coordinates for a given figure size. The result is in the same
% unit as the figure unit. 
% ratio is a 1x2 vector: height:width (default: 3 4)
% margins is a 1x4 vector: left, right, top, bottom
% default: [0.05 0.01 0.05 0.05]
% wmargins is the margins between plots of each axes 1 x 2


switch nargin
    case 6
        wmargins = [0.02 0.02];
    case 5
        margins = [0.05 0.01 0.05 0.05];
        wmargins = [0.02 0.02];
    case 4
        ratio = [3 4];
        margins = [0.05 0.01 0.05 0.05];
        wmargins = [0.02 0.02];
end

% [l_margin, r_margin, t_margin, b_margin]
l_margin = margins(1);
r_margin = margins(2);
t_margin = margins(3);
b_margin = margins(4); 
% percentage of margin between plots of each axes
% pw = 0.08;
% ph = 0.02; % for BM
pw = wmargins(1);
ph = wmargins(2);

x_w = (1-l_margin-r_margin-pw)/nc;
y_h = (1-b_margin-t_margin-ph)/nr;

xspace = (1-nc*x_w-l_margin-r_margin)/(nc-1); % spacing between axes excluding first/last
yspace = (1-nr*y_h-b_margin-t_margin)/(nr-1); % spacing between axes excluding first/last

x_st = [];
x_st(1) = l_margin; 
for x = 2:nc
    x_st(x) = x_st(x-1)+ x_w + xspace;
end

y_st = [];
y_st(1) = b_margin;
for y = 2:nr 
    y_st(y) = y_st(y-1) + y_h + yspace;
end
% y coordinates descend
y_st = fliplr(y_st);

% convert to figure unit
x_st = fw.*x_st;
y_st = fh.*y_st;
x_w = fw*x_w; y_h = fh*y_h;

if ratio(1) > ratio(2) && x_w > y_h % portrait
    x_w = (ratio(2)/ratio(1))*y_h;
elseif ratio(1) > ratio(2) && x_w < y_h
    if y_h > (ratio(1)/ratio(2))*x_w
        y_h = (ratio(2)/ratio(1))*x_w;
    elseif y_h < (ratio(1)/ratio(2))*x_w
        x_w = (ratio(2)/ratio(1))*y_h;
    end
elseif ratio(1) < ratio(2) && x_w > y_h
    if x_w > (ratio(2)/ratio(1))*y_h
        x_w = (ratio(2)/ratio(1))*y_h;
    elseif x_w < (ratio(2)/ratio(1))*y_h
        y_h = (ratio(1)/ratio(2))*x_w;
    end
elseif ratio(1) < ratio(2) && x_w < y_h
    y_h = (ratio(1)/ratio(2))*x_w;
end
return