classdef WordsTrie < handle
    properties
        trie = [];
    end
    
    methods
        function WT = WordsTrie()
            % Create a base node
            node.value = '';
            node.children = [];
            node.is_word = 0;
            WT.trie = node;
            return
        end
        
        function print(WT)
            string = blanks(length(WT.trie)-1);
            for i=1:length(WT.trie)-1
                string(i) = WT.trie(i+1).value;
            end
            disp(string)
        end
        
        function print_words(WT)
            % Build words by children of first node
            str = '';
            ID = 1;
            words = {};
            words = WT.build_words(str, ID, words);
            N = length(words);
            str = '';
            lines = 0;
            for i = 1:N
                str = [str, words{i}, ', '];
                if length(str)>100*(1+lines)
                    lines = lines+1;
                    str = [str(1:end-2),char(10)];
                end
            end
            disp(str(1:end-2))
            disp(['Words in dictionary: ',num2str(N)])
        end
        
        function print_tree(WT)
            N = length(WT.trie);
            tree = zeros(1,N);
            letters = blanks(N);
            is_word = zeros(1,N);
            
            [tree, letters, is_word] = print_tree_step(WT, 1, ...
                    tree, letters, is_word);
                
            treeplot(tree);
            [x,y] = treelayout(tree);
            x = x';
            y = y';
            hold on
            scatter(x(is_word == 1),y(is_word == 1),250,'o','fill', ...
                'MarkerFaceColor',[0 0.8 0], ...
                'MarkerEdgeColor',[0 0 0], 'LineWidth', 2)
            scatter(x(is_word == 0),y(is_word == 0),250,'o','fill', ...
                'MarkerFaceColor',[1 1 1], ...
                'MarkerEdgeColor',[0 0 0], 'LineWidth', 2)
            text(x(:,1), y(:,1), cellstr(letters'), ...
                'FontSize', 12, 'VerticalAlignment','middle', ...
                'HorizontalAlignment','center')
            axis off
            set(gcf,'Color',[1 1 1])
            set(findobj('Type','line'),'Color',[0 0 0],'LineWidth',2)
        end
        
        function [tr, lt, wr] = print_tree_step(WT, ID, tr, lt, wr)
            % Receive a node ID
            children = WT.trie(ID).children;
            
            for i = 1:length(children)
                % Save ID as the parent
                tr(children(i)) = ID;
                lt(children(i)) = WT.trie(children(i)).value;
                wr(children(i)) = WT.trie(children(i)).is_word;
                
                [tr, lt, wr] = print_tree_step(WT, children(i), ...
                    tr, lt, wr);
            end
        end
            
        function words = build_words(WT, str, ID, words)
            children = WT.trie(ID).children;
            N = length(children);
            for i = 1:N
                n_ID = children(i);
                n_str = [str, WT.trie(n_ID).value];
                if WT.trie(n_ID).is_word
                    words = [words, n_str];
                end
                words = build_words(WT, n_str, n_ID, words);
            end
        end
            
            
        function res = add_word(WT, word)
            % See if the word exists
            res = WT.find(word);
            
            if res{1} == 0
                % Word not found in trie, add it as a child of the last
                % related node
                ID = res{2};
                N = length(res{3});
                for i = 1:N
                    new_ID = length(WT.trie)+1;
                    node.value = res{3}(i);
                    node.children = [];
                    node.is_word = 0;
                    WT.trie(new_ID) = node;
                    WT.trie(ID).children = [WT.trie(ID).children, ...
                                                new_ID];
                    ID = new_ID;
                end
                WT.trie(end).is_word = 1;
                disp([word, ' was added successfully']);
                res = 1;
                return
            end
            if res{1} == 1
                % String was already in dictionary but not as a word
                WT.trie(res{2}).is_word = 1;
                disp(['String ''',word,''' set as word']);
                res = 1;
                return
            end
            if res{1} == 2
                % Word was already in dictionary
                disp('Word was already in dictionary');
                res = 0;
                return
            end
        end
            
        function res = find(WT, word)
            % Divide the word into characters
            ID = 1;
            c = 1;
            N = length(word);
            while true
                % Check the c'th level
                res = WT.is_child(word(c), ID);
                if res(1) == 0
                    % The word wasn't found
                    % Return the last found node and remaining word piece
                    if c == 1
                        res = {0, 1, word};
                    else
                        res = {0, ID, word(c:end)};
                    end
                    return
                end
                if res(1) == 1
                    % The letter was found
                    if c == N
                        if WT.trie(res(2)).is_word
                            % The word was found
                            res = {2, word};
                        else
                            % The string was found
                            res = {1, res(2)};
                        end
                        return
                    else
                        % Go deeper
                        ID = res(2);
                        c = c + 1;
                        continue
                    end
                end
            end
        end
        
        function res = is_child(WT, letter, ID)
            % Get children of ID node
            for i=1:length(WT.trie(ID).children)
                % Check if the first letter or all of value is a match
                child_ID = WT.trie(ID).children(i);
                child = WT.trie(child_ID).value;
                if strcmp(letter, child) == 1
                    % The letter is a child of ID
                    res = [1, child_ID];
                    return
                end
            end
            
            res = 0;
        end
    end
end
    