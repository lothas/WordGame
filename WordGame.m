classdef WordGame < handle
    properties
        letters = '';
        rows = 1;
        cols = 1;
        grid = '';
        
        levels = {};
        lvls_filename = 'WG_levels.mat';
        
        sol_words = {};
        sol_len = [];
        sol_moves = [];
        hint = [1,1];
        
        dict = [];
        dict_filename = 'words.mat';
        V = [];     % Vision object
        
        NN = [];    % Neural network object
        NN_filename = 'train_data.mat'; % NN training data file
        NN_X = [];  % Training samples for the NN
        NN_T = [];  % Training targets for the NN
        
        motions = {[-1,-1],[-1,0],[-1,1],...
                    [0,-1],        [0,1],...
                    [1,-1], [1,0], [1,1]};
                
        % Properties for UI
        fig = [];
        uig_axes = [];
        uiw_axes = [];
        im_axes = [];
        grid_handles = [];
        text_handles = [];
        word_handles = [];
        word_text_handles = [];
        word_rows = [];
        ui_word = '';
        ui_pressed = {};
        
        drag = false;
        grid_vis = [];
        ctrl_pressed = false;
        
        gap = 0.02;
        grid_pos = [0, 0, 1, .7];
        row_pos = 0;
        col_pos = 0;
        max_tiles = 15;
        bg_color = [0.18 0.18 0.18];
        color1 = [0.96, 0.96, 0.96];
        color2 = [0.3, 0.28, 0.26];
        color3 = [0.75, 1.00, 0.75];
        color4 = [0.8 0.6 0];
    end
    
    methods
        % Class constructor
        function WG = WordGame(rows, cols, letters)
            % Load dictionary
            WG.load_dict(WG.dict_filename)
            
            % Load levels
            try
                data = load(WG.lvls_filename);
                WG.levels = data.WG_levels;
            catch
                warning('No saved levels found')
            end
            
            % Initialize vision
            try
                WG.V = WGVision('laptop');
            catch
                WG.V = WGVision('office');
            end
            WG.V.capture_cb = @WG.capture_letters;
            
            % Load character recognition NN
            try
                data = load('NN.mat');
                WG.NN = data.net1;
                WG.V.des_size = data.des_size;
            catch
                Warning(['NN not found.',...
                    ' Character recognition unavailable']);
            end
            
            if nargin<3
                if ~isempty(WG.levels)
                    WG.set_game(WG.levels{100,1},WG.levels{100,2});
                else
                    WG.set_game('worbrdain',{'word','brain'});
                end
                return
            end

            WG.set_grid(rows, cols, letters);
        end
        
        %%%%%%%%%%%%%%%%%%% Dictionary functions %%%%%%%%%%%%%%%%%%%
        % %%% FUNCTION: Load Dictionary
        % Loads the dictionary saved in 'filename'.
        function load_dict(WG, filename)
            try
                data = load(filename);
                WG.dict = data.dict;
            catch
                Warning(['''',filename,...
                    ''' not found, new dictionary initialized']);
                WG.dict = WordsTrie();
            end
        end
        
        % %%% FUNCTION: Save Dictionary
        % Saves the dictionary into 'filename'.
        function save_dict(WG)
            dict = WG.dict; %#ok<PROP,NASGU>
            save(WG.dict_filename, 'dict');
        end
        
        % %%% FUNCTION: Add Word
        % Adds the provided word to the dictionary.
        function add_word(WG, word, save_flag)
            if nargin<3
                save_flag = true;
            end
            WG.dict.add_word(word);
            if save_flag
                % Save updated dictionary to file.
                WG.save_dict();
            end
        end
        
        % %%% FUNCTION: Add Words
        % Adds the provided words to the dictionary.
        function add_words(WG, varargin)
            for i = 1:nargin-1
                WG.dict.add_word(varargin{i});
            end
            % Save updated dictionary to file.
            WG.save_dict();
        end
        
        %%%%%%%%%%%%%%%%%%%%%% Game functions %%%%%%%%%%%%%%%%%%%%%%
        % %%% FUNCTION: Set Game
        % Receives an array of chars (grid row by row) and a cell of
        % solutions (words).
        function set_game(WG, letters, words)
            cur_grid = [WG.rows, WG.cols];
            cur_wor = WG.sol_len;
            
            % Set the game grid
            WG.set_grid(letters);
            
            if iscell(words)
                % If "words" provided as solutions:
                WG.sol_words = words;
                WG.sol_len = cellfun(@length,words);
            else
                % If "words" provided as solution length
                WG.sol_len = words;
                % Initialize empty solution cell
                WG.sol_words = cell(1,length(words));
                WG.sol_words(:) = {''};
            end
            WG.sol_moves = cell(1,length(words));
            if length(letters) ~= sum(WG.sol_len)
                error(['ERROR! Word length doesn''t ', ...
                       'match number of letters']);
            end
            
            % Prepare word grid
            m_tiles = WG.max_tiles - 2;
            n_rows = 1;
            WG.word_rows = zeros(length(WG.sol_len),1);
            fit = false;
            while ~fit
                fit = true;
                len = 0;
                row = 1;
                for w = 1:length(WG.sol_len)
                    len = len+WG.sol_len(w);
                    if len <= m_tiles
                        len = len+1; % Add a space
                    else
                        % Word doesn't fit in row, go to next row
                        row = row+1;
                        if row>n_rows
                            % Too many rows, increase max width and
                            % start again
                            m_tiles = m_tiles + 1;
                            n_rows = n_rows + 1;
                            fit = false;
                            break
                        end
                        WG.word_rows(w) = row;
                        len = WG.sol_len(w)+1;
                    end
                    WG.word_rows(w) = row;
                end
            end
                
            % Reset hints
            WG.hint = [1,1];
            
            % Draw the game
            if any([WG.rows, WG.cols]~=cur_grid) || ...
                    length(WG.sol_len)~=length(cur_wor) || ...
                    any(WG.sol_len~=cur_wor)
                close(WG.fig)
            end
            WG.draw_game();
        end
        
        % %%% FUNCTION: Set Grid
        % Receives an array of chars (grid row by row) and dimensions
        % (optional).
        function set_grid(WG, varargin)
            if nargin == 4
                % Letters and dimensions provided
                WG.rows = varargin{1};
                WG.cols = varargin{2};
                WG.letters = varargin{3};
            end
            if nargin == 2
                % Only letters provided, use square grid
                r = sqrt(length(varargin{1}));
                if mod(r,1) ~= 0
                    error('ERROR! Grid provided is not square');
                end
                WG.letters = varargin{1};
                WG.rows = r;
                WG.cols = r;
            end
            
            % Initialize grid with chars
            WG.grid = cell(WG.rows, WG.cols);
            c = 1;
            for i = 1:WG.rows
                for j = 1:WG.cols
                    WG.grid{i,j} = WG.letters(c);
                    c = c + 1;
                end
            end
            
            % Reset word selection on grid
            WG.ui_word = '';
            WG.ui_pressed = {};
        end
        
        % %%% FUNCTION: Update Grid
        % Receives the current grid and selected tiles.
        % Returns the updated grid after removing the selected tiles.
        function out_grid = update_grid(WG, in_grid, pressed)
            out_grid = in_grid;
            % Set the selected tiles as blank
            for i = 1:length(pressed)
                pos = pressed{i};
                out_grid{pos(1), pos(2)} = '';
            end
            
            % Lower letters wherever there's free space (blank tiles)
            for col = 1:WG.cols
                space = cellfun(@isempty,out_grid(:,col));
                % 0 - taken, 1 - free
                % Find first empty space from bottom
                i = find(space == 1, 1, 'last');
                while true
                    if all(space(1:i) == 1)
                        % If all remaining spaces above are free:
                        break;
                    end
                    if space(i) == 1
                        for j = 1:i-1
                            % Move all items in col down by 1
                            out_grid{i-j+1,col} = ...
                                out_grid{i-j,col};
                            space(i-j+1) = space(i-j);
                            % Update space above
                            out_grid{i-j,col} = '';
                            space(i-j) = 1;
                        end
                    else
                        i = i - 1;
                        if i<=1
                            % If top of grid reached:
                            break;
                        end
                    end
                end
            end
        end
        
        % %%% FUNCTION: Print Grid
        % Displays a text version of the game grid.
        function print_grid(WG)
            for i = 1:WG.rows
                string = ['[ ',blanks(2*WG.cols),']'];
                for j = 1:WG.cols
                    string(2*j+1:2*j+2) = [WG.grid{i,j},' '];
                end
                disp(string)
            end
        end

        %%%%%%%%%%%%%%%%%%%%%% Hint functions %%%%%%%%%%%%%%%%%%%%%%
        % %%% FUNCTION: Get Words
        % Calculates possible solutions of length 'des_length' that begin with
        % string 'first'.
        function varargout = get_words(WG, first, des_length)
            words = {};
            moves = {};
            % Build words from grid
            for i = 1:WG.rows
                for j = 1:WG.cols
                    if isempty(WG.grid{i,j}) || ...
                            strcmp(WG.grid{i,j},first(1)) == 0
                        continue
                    end
                    
                    % Starting from each letter on the grid:
                    pos = [i,j];
                    vis = zeros(size(WG.grid));
                    vis(i,j) = 1;
                    str = WG.grid{i,j};
                    
                    % Build words recursively
                    [new_words, n_moves] = ...
                        WG.grid_step(pos, vis, str, des_length, false);
                    
                    if ~isempty(new_words)
                        words = [words, new_words]; %#ok<AGROW>
                        moves = [moves, n_moves]; %#ok<AGROW>
                    end
                end
            end
            
            % Remove repeated words
            unique_words = unique(words);
            
            % Keep only words starting with 'first'
            fit = strfind(unique_words,first);
            no_first = cellfun(@isempty,fit);
            unique_words(no_first) = [];
            fit(no_first) = [];
            first_first = cellfun(@(x)x(1)==1,fit);
            unique_words(~first_first) = [];
            
            % Search for words in Word dictionary
            words_str = '';
            for i = 1:length(unique_words)
                words_str = [words_str,unique_words{i},' ']; %#ok<AGROW>
            end
            status = dictionary(words_str(1:end-3));
            
            unique_words(~status) = [];
            
            switch nargout
                case 0
                    varargout = {};
                    words_str = '';
                    lines = 0;
                    for i = 1:length(unique_words)
                        words_str = [words_str,unique_words{i},' | ']; %#ok<AGROW>
                        if length(words_str)>(lines+1)*75
                            words_str = [words_str(1:end-3), char(10)];
                            lines = lines+1;
                        end
                    end
                    disp(words_str(1:end-3));
                case 1
                    varargout = unique_words;
                case 2
                    varargout = {words, moves};
            end
        end 
            
        % %%% FUNCTION: Get Hint
        % Calculates possible solutions and displays them as hints.
        % More of the word is shown every time.
        function varargout = get_hint(WG, word_len, hint_letters)
            words = {};
            moves = {};
            % Build words from grid
            for i = 1:WG.rows
                for j = 1:WG.cols
                    if isempty(WG.grid{i,j})
                        continue
                    end
                    
                    % Starting from each letter on the grid:
                    pos = [i,j];
                    vis = zeros(size(WG.grid));
                    vis(i,j) = 1;
                    str = WG.grid{i,j};
                    
                    % Build words recursively
                    [new_words, n_moves] = WG.grid_step(pos, vis, str, word_len);
                    
                    if ~isempty(new_words)
                        words = [words, new_words]; %#ok<AGROW>
                        moves = [moves, n_moves]; %#ok<AGROW>
                    end
                end
            end
            
            % Remove repeated words
            unique_words = unique(words);
            
            switch nargout
                case 0
                    varargout = {};
                    words_str = '';
                    for i = 1:length(unique_words)
                        N = min(word_len,hint_letters);
                        % Format the word. Leave empty spaces as required.
                        % For example: p l a _ _
                        str = blanks(2*N);
                        for j = 1:N
                            str(2*j-1:2*j) = [unique_words{i}(j),' '];
                        end
                        str = [str(1:end-1), ...
                            repmat(' _', 1, word_len - hint_letters)];
                        words_str = [words_str,str,' | ']; %#ok<AGROW>
                    end
                    disp(words_str(1:end-3));
                case 1
                    varargout = unique_words;
                case 2
                    varargout = {words, moves};
            end
        end
        
        % %%% FUNCTION: Solve
        % Finds the solution for the grid (if all words exist in dict.)
        function solve(WG)
            disp('Attempting to solve...')
            t_start = tic;
            
            % Save current grid
            temp_grid = WG.grid;
            % Reset if necessary
            if ~all(all(cellfun(@isempty,WG.grid)==0))
                WG.reset_cb(0, 0, 0);
            end
            
            out = WG.solve_step(WG.grid, 1, t_start);
            
            if out.Solved == 1
                to_solve = find(cellfun(@isempty,WG.sol_moves));
                WG.sol_moves(to_solve) = ...
                    out.Moves(1:length(to_solve));
                disp(out.Words)
            else
                WG.get_all_hints();
            end
            
            % Return to original grid
            WG.grid = temp_grid;
            WG.draw_grid();
        end
        
        function out = solve_step(WG, grid, n, t_start)
            % Save current grid
            temp_grid = WG.grid;
            % Apply passed grid
            WG.grid = grid;
                    
            % Initialize output
            out.Solved = 0; % failed to find a solution
            out.Words = {};
            out.Moves = {};
            
            % See if word n already has a solution
            if ~isempty(WG.sol_words{n})
                % Word already solved, was it removed?
                % (assume words are removed in correct order)
                grid_left = sum(sum(~cellfun(@isempty, WG.grid)));
                if grid_left <= (WG.rows*WG.cols)-sum(WG.sol_len(1:n))
                    % Move on to the next word
                    n_words = WG.sol_words(n);
                    n_moves = WG.sol_moves(n);
                else
                    % Check possible moves for the word
                    len = WG.sol_len(n);
                    [n_words, n_moves] = WG.get_hint(len, len);
                    % Remove words that are not the solution
                    ids = find(cellfun(@(x)strcmp(x,WG.sol_words{n}),...
                                       n_words)==0);
                    n_words(ids) = [];
                    n_moves(ids) = [];
                end
            else
                % Look for possible solutions
                len = WG.sol_len(n);
                [n_words, n_moves] = WG.get_hint(len, len);
            end
            
            if ~isempty(n_words)
                if n == length(WG.sol_len)
                    % That's it! We found a solution
                    out.Solved = 1;
                    out.Words = n_words;
                    out.Moves = n_moves;
                else
                    for i = 1:length(n_words)
                        % Check time elased
                        if toc(t_start) > 10
                            disp('Solution algorithm timed out')
                            out.Solved = -1;
                            out.Words = n_words;
                            out.Moves = n_moves;
                            break
                        end
                            
                        % Take word out of grid
                        this_move = mat2cell(n_moves{i}, ...
                            ones(1,size(n_moves{i},1)),2)';
                        n_grid = WG.update_grid(grid, this_move);

                        % Try next step
                        n_out = WG.solve_step(n_grid, n+1, t_start);
                        
                        if n_out.Solved == 1
                            % That's it! We found a solution
                            out.Solved = 1;
                            out.Words = [n_words{i}, n_out.Words];
                            out.Moves = [n_moves{i}, n_out.Moves];
                            break
                        end
                        
                        if n_out.Solved == -1
                            % Solution timed out
                            out.Solved = -1;
                            out.Words = [n_words{i}, n_out.Words];
                            out.Moves = [n_moves{i}, n_out.Moves];
                            break
                        end
                    end
                end
            end
            
            % Return to original grid
            WG.grid = temp_grid;
        end
        
        % %%% FUNCTION: Get All Hints
        % Calls get_hint for each solution length required
        function get_all_hints(WG)
            disp('Getting all hints...')
            for i = 1:length(WG.sol_len)
                if isempty(WG.sol_words{i})
                    WG.get_hint(WG.sol_len(i), WG.sol_len(i))
                end
            end 
        end
        
        % %%% FUNCTION: Grid Step
        % Receives current grid position, the visited nodes on the
        % grid, the current string and the desired word length.
        % Checks all possible moves: if it can move to the new place
        % in the grid, it calls itself recursively until the desired
        % length is reached. If a word is found it is returned,
        % otherwise it returns {}.
        function [words, moves] = grid_step(WG, pos, vis, str, len, word)
            if nargin<6
                % Search for words, not just strings
                word = true;
            end
            words = {};
            moves = {};
            for i = 1:length(WG.motions)
                new_pos = [pos;
                           pos(end,:)+WG.motions{i}];
                % Check if position is valid
                if all([0,0] < new_pos(end,:) & ...
                        new_pos(end,:) <= [WG.rows, WG.cols])
                    % Move is valid
                    if vis(new_pos(end,1),new_pos(end,2)) == 1
                        % Node already visitied
                        continue
                    end
                    if isempty(WG.grid{new_pos(end,1),new_pos(end,2)})
                        % Node is empty
                        continue
                    end
                    
                    % Mark new position as visited
                    n_vis = vis;
                    n_vis(new_pos(end,1),new_pos(end,2)) = 1;
                    % Add letter to word
                    n_str = [str,WG.grid{new_pos(end,1),new_pos(end,2)}];
                    
                    if word
                        % Check if str is a word
                        res = WG.dict.find(n_str);

                        if res{1} == 0
                            % str is not a substr of a word in the dictionary
                            continue
                        end
                    end
                    
                    if length(n_str) == len
                        % Desired length reached
                        if word
                            % Is str a word?
                            if res{1} == 2
                                % str is a word in dict
                                words = [words, res{2}]; %#ok<AGROW>
                                moves = [moves, new_pos]; %#ok<AGROW>
                            end
                        else
                            % Add string as result
                            words = [words, n_str]; %#ok<AGROW>
                            moves = [moves, new_pos]; %#ok<AGROW>
                        end
                    else
                        % Go deeper
                        [new_words, new_moves] = ...
                            WG.grid_step(new_pos, n_vis, n_str, len, word);
                        if ~isempty(word)
                            words = [words, new_words]; %#ok<AGROW>
                            moves = [moves, new_moves]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%% GUI functions %%%%%%%%%%%%%%%%%%%%%%
        % %%% FUNCTION: Draw Game
        % Renders the game board.
        function draw_game(WG)
            if isempty(WG.fig) || ~ishandle(WG.fig)
                % Open the game window and set callbacks
                WG.fig = figure('KeyPressFcn', @WG.key_down_cb, ...
                       'KeyReleaseFcn', @WG.key_up_cb, ...
                       'WindowButtonDownFcn', @WG.mouse_click, ...
                       'WindowButtonUpFcn', @WG.mouse_release, ...
                       'WindowButtonMotionFcn', @WG.mouse_move);
                set(gcf, 'color', WG.bg_color, 'menu', 'no', ...
                        'units', 'normalized', ...
                        'position', [0.37 0.25 0.25 0.64])
            
                % Add a reset button
                uicontrol('style','pushbutton',...
                    'units','normalized',...
                    'position',[WG.gap, WG.gap, 0.4*(1-3*WG.gap), 0.09],...
                    'string','Reset',...
                    'fontsize',20-WG.cols,'fontweight','bold',...
                    'backgroundcolor',[1 0.86 0.9],...
                    'callback',@WG.reset_cb);

                % Add a hint button
                uicontrol('style','pushbutton',...
                    'units','normalized',...
                    'position',[2*WG.gap+0.4*(1-3*WG.gap), WG.gap,...
                                0.6*(1-3*WG.gap), 0.09],...
                    'string','Get Hint',...
                    'fontsize',20-WG.cols,'fontweight','bold',...
                    'backgroundcolor',[0.86 1 0.9],...
                    'callback',@WG.hint_cb);
            else
                figure(WG.fig);
            end
            
            WG.draw_grid();
            WG.draw_words();
        end
        
        % %%% FUNCTION: Draw Grid
        % Renders the game grid.
        function draw_grid(WG)
            if isempty(WG.grid_handles) || ~ishandle(WG.grid_handles(1,1))
                % Draw grid for first time
                % Calculate tile width and height
                button_w = (1-(WG.cols+1)*WG.gap)/WG.cols;
                button_h = (1-(WG.rows+1)*WG.gap)/WG.rows;
                pos = [WG.gap, 1-button_h-WG.gap, ...
                       button_w, button_h];
                
                % Reset handles and visited grid
                WG.grid_vis = zeros(WG.rows, WG.rows);
                WG.grid_handles = zeros(WG.rows, WG.cols);
                WG.text_handles = zeros(WG.rows, WG.cols);
                                    
                % Preare axes for tile render
                axis_pos = WG.grid_pos;
                axis_pos(2) = 1-axis_pos(2)-axis_pos(4);
                WG.uig_axes = axes;
                axis([0 1 0 1]);
                h = gca;
                set(h,'Position',axis_pos);

                % For each tile:
                for i = 1:WG.rows
                    for j = 1:WG.cols
                        % Draw a rectangle
                        WG.grid_handles(i,j) = ...
                            rectangle('Position',pos,...
                            'FaceColor',WG.color1,...
                            'EdgeColor',WG.color1, 'Curvature',0.2);
                        % Add the tile letter
                        WG.text_handles(i,j) = ...
                            text(pos(1)+pos(3)/2, pos(2)+pos(4)/2, ...
                            upper(WG.grid{i,j}), ...
                            'FontName','Corbel', ...
                            'FontSize', 55-4*WG.cols, ...
                            'FontWeight', 'bold', ...
                            'HorizontalAlignment', 'center');
                        % Update horizontal position
                        pos(1) = pos(1) + button_w + WG.gap;
                    end
                    % Update vertical position
                    pos(2) = pos(2) - button_h - WG.gap;
                    % Reset horizontal position
                    pos(1) = WG.grid_pos(1)+WG.gap;
                end
                axis equal
                axis off
                axis manual
                
                % Preare axes for image render
                axis_pos = WG.grid_pos;
                axis_pos(1) = axis_pos(1)+WG.gap;
                axis_pos(3) = axis_pos(3)-2*WG.gap;
                axis_pos(4) = 1.25*axis_pos(4);
                axis_pos(2) = 1-axis_pos(2)-axis_pos(4);
                WG.im_axes = axes;
                axis([0 1 0 1]);
                h = gca;
                set(h,'Position',axis_pos,'Visible','off');
                
                % Initialize grid positions for mouse move detection
                WG.init_grid();
                % Update grid (hides tiles if grid is initialized with
                % empty tiles)
                WG.draw_grid();
            else
                % Update grid
                axes(WG.uig_axes)
                for i = 1:WG.rows
                    for j = 1:WG.cols
                        if ~isempty(WG.grid{i,j})
                            % Tile is not empty, set it as visible
                            set(WG.grid_handles(i,j),'Visible','on')
                            set(WG.text_handles(i,j),'Visible','on')
                            if WG.grid_vis(i,j) == 0
                                % If not visited show with color 1
                                set(WG.grid_handles(i,j), ...
                                'FaceColor',WG.color1,...
                                'EdgeColor',WG.color1);
                                set(WG.text_handles(i,j), ...
                                'color', [0 0 0]);
                            else
                                if WG.grid_vis(i,j) == 1
                                    % If visited show with color 2
                                    set(WG.grid_handles(i,j), ...
                                    'FaceColor',WG.color2,...
                                    'EdgeColor',WG.color2);
                                    set(WG.text_handles(i,j), ...
                                    'color', WG.color4);
                                else
                                    % If highlighted show with color 3
                                    set(WG.grid_handles(i,j), ...
                                    'FaceColor',WG.color3,...
                                    'EdgeColor',WG.color3);
                                    set(WG.text_handles(i,j), ...
                                    'color', [0 0 0]);
                                end
                            end
                            
                            % Update tile string (necessary when tiles are
                            % removed)
                            set(WG.text_handles(i,j), ...
                                'string', upper(WG.grid{i,j}));
                        else
                            % Tile is empty, set it as hidden
                            set(WG.grid_handles(i,j),'Visible','off')
                            set(WG.text_handles(i,j),'Visible','off')
                        end 
                    end
                end
            end
        end
        
        % %%% FUNCTION: Init Grid
        % Initializes grid positions for mouse move detection.
        function init_grid(WG)
            ui_gap = 2*WG.gap;
            button_w = (WG.grid_pos(3)-(WG.cols+1)*ui_gap)/WG.cols;
            button_h = (WG.grid_pos(4)-(WG.rows+1)*ui_gap)/WG.rows;
            WG.row_pos = zeros(1, 2*WG.rows);
            WG.col_pos = zeros(1, 2*WG.cols);
            pos = [WG.grid_pos(1)+ui_gap, ...
                   1-WG.grid_pos(2)-button_h-ui_gap, ...
                   button_w, button_h];
            
            for i = 1:WG.rows
                WG.row_pos(2*i-1:2*i) = [pos(2)+pos(4), pos(2)];
                for j = 1:WG.cols
                    WG.col_pos(2*j-1:2*j) = [pos(1), pos(1)+pos(3)];
                    pos(1) = pos(1) + button_w + ui_gap;
                end
                pos(2) = pos(2) - button_h - ui_gap;
                pos(1) = WG.grid_pos(1)+ui_gap;
            end
        end
        
        % %%% FUNCTION: Draw Words
        % Draws the solved/unsolved words.
        function draw_words(WG)
            if isempty(WG.word_handles) || ~ishandle(WG.word_handles(1,1))
                % Draw words for first time
                
                % Preare axes for tile render
                axis_pos = WG.grid_pos;
                axis_pos(2) = 0.12;
                axis_pos(4) = 1-axis_pos(2)-axis_pos(4)-WG.grid_pos(2);
                WG.uiw_axes = axes;
                axis([0 1 0 1]);
                h = gca;
                axis_pos(1) = axis_pos(1)+0.06;
                axis_pos(3) = axis_pos(3)-0.12;
                set(h,'Position',axis_pos);
                AR = axis_pos(4)/axis_pos(3);
                
                % Calculate tile width and height
                tgap = 0.009;
                n_rows = max(WG.word_rows);
                h_gap = (n_rows-1)*2*tgap;
                n_tiles = zeros(n_rows,1);
                for i = 1:n_rows
                    n_tiles(i) = sum(WG.sol_len(WG.word_rows==i)) + ...
                        length(WG.sol_len(WG.word_rows==i)) - 1;
                end
                button_w = (1-(max(n_tiles)-1)*tgap-2*h_gap)/(max(n_tiles));
%                 button_w = (1 + tgap - 2*tgap)/max(n_tiles) - tgap;
                button_h = (button_w-tgap)/AR;
                d_center = (1-n_rows*button_h-(n_rows-1)*tgap/AR)/2;
                pos = [h_gap, 1-button_h-d_center, ...
                       button_w, button_h];
                
                wh = 1;
                WG.word_handles = [];
                WG.word_text_handles = [];
                for i = 1:n_rows
                    pos(1) = h_gap + (button_w+tgap) * ...
                                      (max(n_tiles)-n_tiles(i))/2;
                    words = find(WG.word_rows == i);
                    for w = 1:length(words)
                        if isempty(WG.sol_moves{words(w)})
                            word_s = blanks(WG.sol_len(words(w)));
                        else
                            word_s = WG.sol_words{words(w)};
                        end
                        
                        for l = 1:WG.sol_len(words(w))
                            % Add a word tile
                            if isletter(word_s(l))
                                WG.word_handles(wh) = ...
                                rectangle('Position',pos,...
                                'FaceColor',WG.color1,...
                                'EdgeColor',WG.color1, 'Curvature',0.5);
                            else
                                WG.word_handles(wh) = ...
                                rectangle('Position',pos,...
                                'FaceColor',WG.bg_color,...
                                'EdgeColor',WG.color4, 'Curvature',0.5);
                            end
                            
                            % Add the tile letter
                            WG.word_text_handles(wh) = ...
                                text(pos(1)+pos(3)/2, pos(2)+pos(4)/2, ...
                                upper(word_s(l)), ...
                                'fontsize', 14-n_rows, ...
                                'fontweight', 'bold', ...
                                'HorizontalAlignment', 'center');
                            
                            % Update horizontal position
                            pos(1) = pos(1) + button_w + tgap;
                            wh = wh + 1;
                        end
                        
                        % Add a space between words
                        pos(1) = pos(1) + button_w + tgap;
                    end
                    
                    % Update vertical position
                    pos(2) = pos(2) - button_h - tgap/AR;
                end
                
                axis off
                axis manual
            else
                % Update grid
                axes(WG.uiw_axes)
                
                n_rows = max(WG.word_rows);
                wh = 1;
                for i = 1:n_rows
                    words = find(WG.word_rows == i);
                    for w = 1:length(words)
                        if isempty(WG.sol_moves{words(w)})
                            word_s = blanks(WG.sol_len(words(w)));
                        else
                            word_s = WG.sol_words{words(w)};
                            if isempty(word_s)
                                word_s = blanks(WG.sol_len(words(w)));
                            end
                        end
                        
                        for l = 1:WG.sol_len(words(w))
                            if isletter(word_s(l))
                                set(WG.word_handles(wh), ...
                                'FaceColor',WG.color1,...
                                'EdgeColor',WG.color1);
                            else
                                set(WG.word_handles(wh), ...
                                'FaceColor',WG.bg_color,...
                                'EdgeColor',WG.color4);
                            end
                            set(WG.word_text_handles(wh), ...
                            'String',upper(word_s(l)));
                        
                            wh = wh + 1;
                        end
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%% GUI interaction functions %%%%%%%%%%%%%%%%%
        % %%% FUNCTION: Mouse Click
        % Callback for mouse click. Starts tile selection.
        function mouse_click(WG, src, data)
            if strcmp(get(WG.fig,'selectiontype'),'alt')
                % Right button was clicked
                mouse_xy = get(WG.fig,'CurrentPoint');
                
                ax_pos = get(WG.uiw_axes, 'Position');
                if mouse_xy(2) > ax_pos(2) && ...
                        mouse_xy(2) < ax_pos(2) + ax_pos(4)
                    % Right click on word solution area
                    res = WG.myinput(...
                        'Enter new solution lengths separated by commas', ...
                        'Change solution lengths');
                    if isempty(res)
                        return
                    end
                    
                    let = WG.letters;
                    res = res(regexp(res,'[0-9,]'));
                    % accept only numbers separated by commas
                    len = str2num(res); %#ok<ST2NM>
                    
                    WG.set_game(let, len);
                    return
                end 
                
                % Check if mouse position coincides with a tile
                grid_xy = WG.get_button_by_xy(mouse_xy);

                if ~isempty(grid_xy)
                    % Mouse is over a letter
                    res = WG.myinput(...
                        'Enter a new character', ...
                        'Change tile');
                    if isempty(res)
                        return
                    end
                    ind = sub2ind(size(WG.grid), grid_xy(2), grid_xy(1));
                    WG.grid{grid_xy(1), grid_xy(2)} = lower(res(1));
                    WG.letters(ind) = lower(res(1));
                    
                    % NN output was corrected, add it to the training data
                    WG.NN_X = [WG.NN_X, WG.V.letters(:,ind)];
                    T = zeros(26,1); T(uint8(lower(res(1)))-96) = 1;
                    WG.NN_T = [WG.NN_T, T];
                    
                    WG.draw_game();
                    return
                end
            end
            
            WG.drag = true;
            mouse_move(WG, src, data);
        end
        
        % %%% FUNCTION: Mouse Release
        % Callback for mouse release. Ends tile selection and checks if
        % selection is a word.
        function mouse_release(WG, ~, ~)
            WG.drag = false;
            if ~isempty(WG.ui_word)
                WG.process_word();
            end
        end
        
        function process_word(WG)
            is_word = false;
            
            res = WG.dict.find(WG.ui_word);
            if res{1} == 2
                is_word = true;
                % Check if a word of that length is needed
                N = length(WG.ui_word);

                % Words of length N
                pos_len = WG.sol_len == N;
                % Words which have been found
                pos_mov = ~cellfun(@isempty,WG.sol_moves);
                pos_wor = ~cellfun(@isempty,WG.sol_words);
                % Is the word a found/provided solution?
                pos_sol = strcmp(WG.ui_word,WG.sol_words);

                pos = find(pos_sol == 1, 1, 'first');
                if ~isempty(pos)
                    % Word is already marked as a solution
                    % Has it been taken out of the grid?
                    if ~pos_mov(pos)
                        % Mark it as solution
                        WG.sol_moves{pos} = cell2mat(WG.ui_pressed');
                        WG.sol_words{pos} = WG.ui_word;
                    else
                        % Word already taken out of grid
                        is_word = false;
                    end
                else
                    % Word doesn't appear as solution yet
                    pos = find(pos_len & ~pos_mov & ~pos_wor, 1, 'first');
                    if ~isempty(pos)
                        % A word of that length is needed

                        % Mark it as solution
                        WG.sol_moves{pos} = cell2mat(WG.ui_pressed');
                        WG.sol_words{pos} = WG.ui_word;
                    else
                        % A word of that length isn't needed
                        is_word = false;
                    end
                end
            end

            if is_word
                % Word found, take it out of grid
                WG.grid = update_grid(WG, WG.grid, WG.ui_pressed);
                WG.draw_grid();
                WG.draw_words();
            end 
            
            % Return tiles to original color
            WG.grid_vis = 0*WG.grid_vis;
            WG.draw_grid();
            
            WG.ui_word = '';
            WG.ui_pressed = {};
        end

        % %%% FUNCTION: Mouse move
        % Callback for mouse motion. Selects tiles when mouse is dragged
        % over them while pressed.
        function mouse_move(WG, ~, ~)
            if ~WG.drag
                return
            else
                mouse_xy = get(WG.fig,'CurrentPoint');
                
                % Check if mouse position coincides with a tile
                grid_xy = WG.get_button_by_xy(mouse_xy);
                
                new = false;
                if ~isempty(grid_xy)
                    % Mouse is over a letter
                    if isempty(WG.ui_pressed)
                        % This is the first letter selected
                        new = true;
                    else
                        if WG.grid_vis(grid_xy(1), grid_xy(2)) == 0
                            % We haven't visited this letter
                            new = true;
                        end
                    end
                end
                
                if new
                    % Add the letter to the selected word
                    WG.ui_word = [WG.ui_word, ...
                        WG.grid{grid_xy(1), grid_xy(2)}];
                    % Add the tile to selected tiles
                    WG.ui_pressed = [WG.ui_pressed, grid_xy];
                    % Mark tile as visitied
                    WG.grid_vis(grid_xy(1), grid_xy(2)) = 1;
                    % Update grid display
                    WG.draw_grid();
                end
            end
        end
        
        % %%% FUNCTION: Get Button By X Y
        % Checks if the mouse position coincides with a tile.
        function pos = get_button_by_xy(WG, mouse_xy)
            col = find(mouse_xy(1)<WG.col_pos,1,'first');
            row = find(mouse_xy(2)>WG.row_pos,1,'first');
            if isempty(col) || isempty(row) || any(mod([row,col],2))
                % Mouse is in a gap
                pos = [];
            else
                % Mouse is over a tile
                if ~isempty(WG.grid{row/2, col/2})
                    pos = [row/2, col/2];
                else
                    pos = [];
                end
            end
        end
        
        function key_down_cb(WG, ~, ~, ~)
            k = get (gcf, 'CurrentKey');
            
            if strcmp(k,'control')
                WG.ctrl_pressed = true;
                return
            end
            
            if strcmp(k,'space')
                imh = findobj(WG.im_axes,'Type','Image');
                if ~isempty(imh)
                    set(imh,'Visible','on');
                    uistack(WG.im_axes,'top')
                end
%                 WG.draw_game();
                return
            end
            
            % Show words of certain length that begin with certain string
            if k(1) == 'f' && WG.ctrl_pressed
                first = WG.myinput(...
                    'Enter the first letter or letters', ...
                    'Starts with...');
                if isempty(first) || any(~isletter(first))
                    return
                end
                
                des_length = WG.myinput(...
                    'Enter the desired word length', ...
                    'Length');
                if isempty(des_length) || isnan(str2double(des_length))
                    return
                end
                
                WG.get_words(first, str2double(des_length));
            end
            
            % Show all hints
            if k(1) == 'x' && WG.ctrl_pressed
                WG.get_all_hints();
                set(gcf,'CurrentCharacter','a');
                return
            end
            
            % Attempt to solve
            if k(1) == 's' && WG.ctrl_pressed
                WG.solve();
                set(gcf,'CurrentCharacter','a');
                return
            end
            
            % Reset
            if k(1) == 'r' && WG.ctrl_pressed
                WG.reset_cb(0, 0, 0);
                set(gcf,'CurrentCharacter','a');
                return
            end
            
            % Open "Add words" dialog
            if k(1) == 'a' && WG.ctrl_pressed
                res = WG.myinput(...
                    'Enter a word (or words separated by commas)', ...
                    'Add word(s)');
                if isempty(res)
                    return
                end
                
                % Process input
                res(regexp(res,'[ ]'))=[]; % clean spaces
                words = strsplit(res,','); % split by commas
                for i = 1:length(words)
                    w = lower(words{i});

                    % Check that word is valid (has only letters a-z)
                    if strcmp(w,w(regexp(w,'[a-z]')))
                        WG.add_word(w);
                    end
                end
                
                return
            end
            
            % Open "set game" dialog
            if k(1) == 'n' && WG.ctrl_pressed
                if isempty(WG.NN)
                    res = WG.myinput(...
                        'Enter new grid letters row by row', ...
                        'Set new game grid');
                    if isempty(res)
                        return
                    end
                
                    % Process input
                    res = lower(res);
                    let = res(regexp(res,'[a-z]')); % accept only letters
                    
                    res = WG.myinput(...
                        ['Enter the game solution words (or their', ...
                         'length) separated by commas'], ...
                         'Set new game words');
                    if isempty(res)
                        return
                    end

                    if sum(isletter(res))>0.5*length(res)
                        % Solutions inserted as words
                        res = lower(res);
                        res = res(regexp(res,'[a-z,]'));
                        % accept only words separated by commas
                        words = strsplit(res,','); % split by commas
                    else
                        % Solutions inserted as length
                        res = res(regexp(res,'[0-9,]'));
                        % accept only numbers separated by commas
                        words = str2num(res); %#ok<ST2NM>
                    end
                    WG.set_game(let,words);
                else
                    % Save NN training data before setting a new game
                    if ~isempty(WG.NN_X)
                        try
                            data = load(WG.NN_filename);
                            NN_samples = [data.NN_samples, WG.NN_X]; %#ok<NASGU>
                            NN_targets = [data.NN_targets, WG.NN_T]; %#ok<NASGU>
                        catch
                            NN_samples = WG.NN_X; %#ok<NASGU>
                            NN_targets = WG.NN_T; %#ok<NASGU>
                        end
                        save(WG.NN_filename,'NN_samples','NN_targets');
                        WG.NN_X = [];
                        WG.NN_T = [];
                    end
                    
                    % Save game if solved
                    if ~any(cellfun(@isempty,WG.sol_moves))
                        WG_levels = [WG.levels;
                                     WG.letters, {WG.sol_words}];
                        save(WG.lvls_filename,'WG_levels');
                        WG.levels = WG_levels;
                        WG.sol_moves = cell(1,length(WG.sol_words));
                    end
                    
                    WG.V.start_capture();
                    waitfor(WG.V.fh);
                    
                    axes(WG.im_axes)
                    h = image(repmat(WG.V.I_store{1},[1,1,3])/255);
                    set(h, 'Visible', 'off');
                    set(WG.im_axes,'Visible','off')
%                     alpha = 0.5*ones(size(WG.V.I_store{1}));
%                     set(h, 'AlphaData', alpha, 'Visible', 'off');
                end
                
                return
            end
            
            if k(1) == '1'
                tiles_out = sum(sum(cellfun(@isempty,WG.grid)==1));
                next = 1; sum_len = 0;
                while tiles_out > sum_len
                    sum_len = sum_len + WG.sol_len(next);
                    next = next+1;
                end
                
                if ~isempty(next)
                    % Highlight solution moves
                    moves = WG.sol_moves{next};
                    for m = 1:size(moves,1)
                        WG.grid_vis(moves(m,1),moves(m,2)) = 2;
                    end
                    WG.draw_grid();
                end
                return
            end
            
            if k(1) == '2'
                tiles_out = sum(sum(cellfun(@isempty,WG.grid)==1));
                next = 1; sum_len = 0;
                while tiles_out > sum_len
                    sum_len = sum_len + WG.sol_len(next);
                    next = next+1;
                end
                
                if ~isempty(next) && next<=length(WG.sol_len) ...
                        && ~isempty(WG.sol_moves{next})
                    % Add word as a solution
                    moves = WG.sol_moves{next};
                    word = blanks(WG.sol_len(next));
                    for m = 1:size(moves,1)
                        word(m) = WG.grid{moves(m,1),moves(m,2)};
                    end
                    WG.sol_words{next} = word;
                    
                    % Take it out of grid
                    this_move = mat2cell(moves, ...
                        ones(1,size(WG.sol_moves{next},1)),2)';
                    WG.grid = WG.update_grid(WG.grid, this_move);
                    WG.draw_game();
                end
                return
            end
            
            % Highlight selected character
            for i = 1:WG.rows
                for j = 1:WG.cols
                    if WG.grid{i,j} == k
                        WG.grid_vis(i,j) = 2;
                    end
                end
            end
            WG.draw_grid();
        end
        
        function key_up_cb(WG, ~, ~, ~)
            k = get (gcf, 'CurrentKey');
            
            if strcmp(k,'control')
                WG.ctrl_pressed = false;
                return
            end
            
            if strcmp(k,'space')
                imh = findobj(WG.im_axes,'Type','Image');
                if ~isempty(imh)
                    set(imh,'Visible','off');
                    uistack(WG.im_axes,'bottom')
                end
%                 WG.draw_game();
                return
            end
            
            if length(k)==1 && (isletter(k(1)) || str2double(k(1))==1)
                % Return characters to normal
                for i = 1:WG.rows
                    for j = 1:WG.cols
                        if any(cellfun(@(x)all(x==[i,j]),WG.ui_pressed))
                            WG.grid_vis(i,j) = 1;
                        else
                            WG.grid_vis(i,j) = 0;
                        end
                    end
                end
                WG.draw_grid();
                drawnow;
            end
        end
        
        function capture_letters(WG, WGVh)
            X = WG.NN(WGVh.letters);    
            R = size(X,2);
            res = blanks(R);
            for n = 1:R
                res(n) = char(96+find(X(:,n)==max(X(:,n))));
            end
            close(WG.V.fh);
            WG.set_game(res,WGVh.word_len);
            
            % Show # of samples of letters captured
            data = load(WG.NN_filename);
            ltrs = unique(res);
            u = blanks(4*length(ltrs));
            d = u;
            ns = 3; % spaces for each char display
            for i = 1:length(ltrs)
                c = uint8(ltrs(i))-96;
                Ns = sum(data.NN_targets(c,:) == 1);
                if Ns <= 35
                    u(1+ns*(c-1):ns*c) = ...
                        sprintf([blanks(ns-1),'%c'], char(96+c));
                    d(1+ns*(c-1):ns*c) = ...
                        sprintf(['%',num2str(ns),'d'], Ns);
                end
            end
            % Show samples per letter
            disp(u); disp(d);
            
%             % Show image and result side by side
%             gcf_pos = get(WG.fig, 'Position');
%             gcf_units = get(WG.fig, 'Units');
%             figure('menu', 'no', 'units', gcf_units, ...
%                 'Position',[gcf_pos(1)+1.1*gcf_pos(3), ...
%                             gcf_pos(2)+0.2*gcf_pos(4), ...
%                             gcf_pos(3), gcf_pos(4)*0.8]);
%             imshow(double(WG.V.I_store{1})/255);
        end
        
        function reset_cb(WG, ~, ~, ~)
            if all(all(cellfun(@isempty,WG.grid)==0))
                % Second reset, clear game and solutions
                WG.set_game(WG.letters, WG.sol_len);
            else
                % Clear grid and moves (but not solution words)
                WG.set_grid(WG.letters);
                WG.sol_moves = cell(1,length(WG.sol_words));
                WG.draw_game();
            end
        end
        
        function hint_cb(WG, ~, ~, ~)
            to_solve = find(cellfun(@isempty,WG.sol_words),1,'first');
            if WG.hint(1) < to_solve
                WG.hint = [to_solve, 1];
            end
            
            WG.get_hint(WG.sol_len(WG.hint(1)),WG.hint(2));
            WG.hint(2) = WG.hint(2)+1;
            if WG.hint(2) > WG.sol_len(WG.hint(1)) ...
                    && WG.hint(1) < length(WG.sol_len)
                WG.hint(1) = WG.hint(1) + 1;
                WG.hint(2) = 1;
            end
        end
        
        function train_net(WG)
            % Load training data
            data = load(WG.NN_filename);
            
            % Get the number of samples for each letter
            u = blanks(4*26);
            d = u;
            ns = 3; % spaces for each char display
            for c = 1:26
                samples = find(data.NN_targets(c,:) == 1);
                u(1+ns*(c-1):ns*c) = sprintf([blanks(ns-1),'%c'], ...
                    char(96+c));
                d(1+ns*(c-1):ns*c) = sprintf(['%',num2str(ns),'d'], ...
                    length(samples));
            end
            % Show samples per letter
            disp(u); disp(d);

            net1 = train(WG.NN, data.NN_samples, data.NN_targets, nnMATLAB);
            des_size = WG.V.des_size; %#ok<NASGU>
            save('NN.mat','net1','des_size');
            WG.NN = net1;
        end
        
        function res = myinput(WG, Prompt, Title) %#ok<INUSL>
            % Input dialog. User may hit return after entering info.

            if nargin<2
                Prompt = 'Enter a word (or words separated by commas)';
            end
            if nargin<3
                Title = Prompt;
            end
            
            width = 480;
            pgap = 10;
            lines = ceil(length(Prompt)/44);

            res = []; % In case the user closes the GUI.
            dlg.fh = figure('units','pixels',...
                          'position',[500 500 width 90+30*lines],...
                          'menubar','none',...
                          'numbertitle','off',...
                          'name',Title,...
                          'resize','off');
            dlg.ed = uicontrol('style','text',...
                             'units','pix',...
                            'position',[pgap 80 width-2*pgap 30*lines],...
                            'string',Prompt, ...
                            'FontSize',16);
            dlg.ed = uicontrol('style','edit',...
                             'units','pix',...
                            'position',[pgap 50 width-2*pgap 30],...
                            'string','', ...
                            'HorizontalAlignment','left');
            dlg.pb = uicontrol('style','pushbutton',...
                             'units','pix',...
                            'position',[width/2-50 10 100 30],...
                            'string','OK', ...
                            'FontSize',16,...
                            'callback',{@pb_call});
            set(dlg.ed,'call',@ed_call)
            uicontrol(dlg.ed) % Make the editbox active.
            uiwait(dlg.fh) % Prevent all other processes from starting until closed.

            function [] = pb_call(varargin)
                res = get(dlg.ed,'string');
                close(dlg.fh); % Closes the GUI, allows the new R to be returned.
            end

            function [] = ed_call(varargin)
                k = get(dlg.fh,'CurrentKey');
                if strcmp(k,'return')
                    uicontrol(dlg.pb)
                    drawnow
                    res = get(dlg.ed,'string');
                    close(gcbf)
                end
            end
        end
    end
end