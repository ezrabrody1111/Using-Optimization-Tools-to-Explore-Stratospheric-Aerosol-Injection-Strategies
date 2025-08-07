function [contours, color_map] = custom_color_map(min,max,contour_levels,colors)
    contours = min:((max-min)/contour_levels):max;
    n_levels = length(contours)-1;
    color_map = brewermap(n_levels,colors);
end