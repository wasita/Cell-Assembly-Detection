function prettifyFig(size, ~)

set(gcf, 'Color', [1 1 1]) % white background
set(gca, 'TickDir', 'out', 'Box', 'off') % flip the y-axis ticks
set(gca, 'FontSize', 20)

if nargin >= 1
    switch size
        case 'sq' % for square (width ~ = length)
            set(gcf,'units','centimeters','position',[0 0 40 40]);
        case 'hb' % for hamburger aka vertical (width > length)
            set(gcf,'units','centimeters','position',[0 0 40 20]);
        case 'hd' % for hotdog aka horizontal (length > width)
            set(gcf,'units','centimeters','position',[0 0 20 40]);
    end
    
    % means there's a colorbar, so want to flip axes there too
    if nargin == 2
        h = colorbar;
        set(h, 'TickDirection', 'out')
    end
    
else
    % default to square
    set(gcf,'units','centimeters','position',[0 0 40 40]);
end

end
