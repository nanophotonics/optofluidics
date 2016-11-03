function x = dialog_table(row_titles, column_titles, default_values)
% function x = dialog_table(default_values)
% X = dialog_table(R,C) allows interactively create arrays of given size instead of typing them in command line. 
% R - number of rows
% C - number of columns

    arr_size = size(default_values);
    titles = cell(size(default_values)+1);
    titles(1,2:end) = column_titles;
    titles(2:end,1)= row_titles;
%     titles{1,1} = '.';

    border_size = 10; % border size
    box_width = 100; % width of the input boxes
    box_height = 20; % height of the input boxes

    fig_dim = [(arr_size(1)+2)*(box_height+border_size)+border_size, ...
               (arr_size(2)+1)*(box_width+border_size)+border_size]; % [y,x]
    arr_fig = figure('unit','pixels',...
                     'NumberTitle','off',...
                     'Menubar','none',...
                     'resize','on',...
                     'position', [300 300 fig_dim(2) fig_dim(1)]);
    ok_but = uicontrol('Style','pushbutton',...
                       'unit','pixels',...
                       'String','Okay',...
                       'position',[border_size border_size fig_dim(2)-2*border_size box_height],...
                       'callback','uiresume',...
                       'tag','ok');
    % cancel_but = uicontrol('Style','pushbutton','unit','pixels','String','Cancel','position',[35 15 40 20],...
    %     'callback','uiresume','tag','cancel');
    

    posx = border_size;
    posxb = posx;
    posy = fig_dim(1)-(box_height+border_size);

    for i = 0:arr_size(1) % rows
        for j = 0:arr_size(2) % colums
            if or(i == 0, j == 0)
                uicontrol('Style','text',...
                          'unit','pixels',...
                          'String', 'row title',...
                          'Position',[posx posy box_width box_height],...
                          'String', titles{i+1,j+1});
            else
                a(i,j) = uicontrol('Style','edit',...
                                   'unit','pixels',...
                                   'String', num2str(default_values(i,j)),...
                                   'position', [posx posy box_width box_height]);
            end
            posx = posx+box_width+border_size;
        end
        posx = posxb;
        posy = posy-box_height-border_size;
    end
    
    

    uiwait;
    but = gco;

    if strcmp(get(but,'tag'),'ok')
        [r,c] = size(a);
        for i = 1:r
            for j = 1:c
                x(i,j) = str2num(get(a(i,j),'string'));
            end
        end
        close;    
    else 
        close;
        return;
end
    
