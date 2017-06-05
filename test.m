function [] = test()

val = fzero(@(x) root(x), [1.1 100])

function [tgap] = root(ar)
    [~, yc, ~, yr] = first_pts(0, ar, 0.3, 1.4);
    
    scale = 0.1622016484252/yc;
    
    yc = yc*scale;
    yr = yr*scale;
    tgap = yc-yr;
    tgap = tgap-0.0196850393700787;
end

end