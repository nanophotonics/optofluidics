function [ colour_RGB ] = colour_gradient(x, x_max, type)

if strcmp(type, 'yellow') == 1
    colour_RGB = [1, 1-x/x_max ,0]; % yellow to red
elseif strcmp(type, 'faded') == 1
    colour_RGB = [1, (1-x/x_max + 0.5)/1.5 , (1-x/x_max)*0.8 ]; % yellow to red - faded
elseif strcmp(type, 'temperature') == 1
    if x/x_max < 0.25
        colour_RGB = [0, x/x_max/0.25, 1]; % temperature
    end
    if x/x_max >= 0.25
        colour_RGB = [0, 1, 1-(x/x_max-0.25)/(0.5-0.25)]; % temperature
    end
    if x/x_max >= 0.5
        colour_RGB = [(x/x_max-0.5)/(0.75-0.5), 1, 0]; % temperature
    end
    if x/x_max >= 0.75
        colour_RGB = [1, 1-(x/x_max-0.75)/(1-0.75), 0]; % temperature
    end
elseif strcmp(type, 'red') == 1
    colour_RGB = [1-x/x_max, 0, 0]; % red
elseif strcmp(type, 'green') == 1
    colour_RGB = [1-x/x_max, 1-x/x_max, 0]; % olive green
elseif strcmp(type, 'aqua') == 1
    colour_RGB = [0, 1-x/x_max, 1-x/x_max]; % aqua
elseif strcmp(type, 'blue') == 1
    colour_RGB = [1-x/x_max, 1-x/x_max, 1]; % blue
elseif strcmp(type, 'purple') == 1
    colour_RGB = [1-x/x_max, 0, 1]; % purple
elseif strcmp(type, 'gray') == 1
    colour_RGB = [1-x/x_max, 1-x/x_max, 1-x/x_max]; % grayscale
end

end

