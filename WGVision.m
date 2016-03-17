classdef WGVision < handle
    %WGVISION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Letters captured
        letters = [];
        word_len = [];
        des_size = [9, 9];  % letter size for NN
        
        % Which hardware are we using?
        laptop = true;
        go_cap = false;
        
        % Video input
        vid = []; % Video hardware
        res = []; % Resolution
        
        % Preview window
        fh = []; % figure handle
        ih = []; % image object handle
        ph = []; % plot handle
        % Base plot data
        pxlim = [];
        pylim = [];
        
        % Image storage
        n_store = 3;
        I_store = {};
        I_n = 1;
        grid_stable = false;
        
        % External function handles
        capture_cb = [];
        done_cb = [];
    end
    
    methods
        function WGV = WGVision(type)
            switch lower(type)
                case {'laptop'}
                    WGV.vid = videoinput('winvideo');
                    WGV.res = get(WGV.vid, 'VideoResolution');
                case {'office'}
                    WGV.laptop = false;

                    WGV.vid = webcam(1);
                    set(WGV.vid, 'Resolution', ...
                        WGV.vid.AvailableResolutions{end-2});
                    WGV.res = cellfun(@str2num, ...
                        strsplit(WGV.vid.AvailableResolutions{end-2},'x'));
            end
            
            % Set plot lims
            WGV.pxlim = [floor(WGV.res(1)*0.75) floor(WGV.res(1)*0.95)];
            WGV.pylim = [floor(WGV.res(2)*0.95) floor(WGV.res(2)*0.75)];
            
            % Prepare image storage
            WGV.I_store = cell(1,WGV.n_store);
        end
        
        function start_capture(WGV)
            %Start preview figure
            WGV.fh = figure('CloseRequestFcn', @WGV.done_capture, ...
                'KeyPressFcn', @WGV.key_press_cb);
            WGV.ih = imshow(zeros(WGV.res(2),WGV.res(1)));
            set(gca,'CLim',[0 255]);
            
            % Prepare plot
            hold on
            WGV.ph = plot(0,0,'LineWidth',3);

            if WGV.laptop
                previewImageObject = preview(WGV.vid, WGV.ih);
                setappdata(previewImageObject, ...
                    'UpdatePreviewWindowFcn', @WGV.update_cb);

                %Set image to be a mirror image view
                Parent = get(previewImageObject,'Parent');
                set(Parent, 'XDir', 'reverse');
            else
                WGV.go_cap = true;
                while WGV.go_cap
                    % Acquire a single image.
                    rgbImage = snapshot(WGV.vid);

                    % Perform image processing
                    WGV.pre_process(rgbImage);
                end
                delete(WGV.fh)
            end
        end
        
        % Update callback for laptop's preview window
        function update_cb(WGV, ~, event, ~)
            % Acquire a single image.
            % image=get(previewImageObject, 'CData');
            rgbImage = event.Data;

            % Perform image processing
            try
                WGV.pre_process(rgbImage);
            catch ME
                disp(['ERROR! ',ME.message])
                disp(ME.cause)
                for i=1:length(ME.stack)
                    disp(['In file: ',ME.stack(i).file,...
                        ', name: ',ME.stack(i).name,...
                        ', line: ',num2str(ME.stack(i).line)])
                end
                delete(gcf)
                rethrow(ME)
            end
        end
        
        function pre_process(WGV, I)            
            % Convert RGB to grayscale.
            grayImage = rgb2gray(I);

            % Detect MSER regions.
            [mserRegions] = detectMSERFeatures(grayImage, ...
                'MaxAreaVariation',0.15, ...
                'ThresholdDelta',20);

            if length(mserRegions)<4
                set(WGV.ih,'CData',grayImage);
                return;
            end
            
            % Calculate region properties
            mserStats = WGV.GetStats(mserRegions, size(grayImage));

            % Take out regions that are not tiles (by comparison)
            
            % Remove regions with very small or very large area
            area = [mserStats.Area];
            ar_mean = mean(area); ar_std = std(area);
            tol = 3;
            area_in = area<ar_mean+tol*ar_std & ...
                area>max(ar_mean-tol*ar_std,90);
            % Keep regions with orientation close to k*90 degrees
%             ori_in = abs(abs([mserStats.Orientation]/90)-0.5)>0.25;
            % Keep regions of low eccentricity
            ecc_in = [mserStats.Eccentricity]<0.6;
            % Keep regions of high solidity
            sol_in = [mserStats.Solidity]>0.4;
            % Keep regions with aspect ratio close to 1
            bb = vertcat(mserStats.BoundingBox);
            asp = double(bb(:,3)')./double(bb(:,4)');
            asp_in = asp>0.7;
            
            reg_in = area_in & ecc_in & sol_in & asp_in;
            mserStats(~reg_in) = [];
            
            if length(mserStats)<4
                set(WGV.ih,'CData',grayImage);
                return;
            end

            % Find bounding box for the grid
            bb = vertcat(mserStats.BoundingBox);
            bb(:,3) = bb(:,1) + bb(:,3);
            bb(:,4) = bb(:,2) + bb(:,4);
            min_bb = min(bb);
            max_bb = max(bb);
            grid_bb = round([min_bb(2) min_bb(1) max_bb(4) max_bb(3)]);

            % Crop image
            try
                Ic = double(grayImage(grid_bb(1):grid_bb(3), grid_bb(2):grid_bb(4)));
                % Get area directly below grid
                min_y = grid_bb(3)+(grid_bb(3)-grid_bb(1))/50;
                max_y = min(grid_bb(3)+(grid_bb(3)-grid_bb(1))/3.3,size(grayImage,1));
                Ic2 = double(grayImage(round(min_y):round(max_y), grid_bb(2):grid_bb(4)));
            catch
                set(WGV.ih,'CData',grayImage);
                set(WGV.ph, 'XData', 0, 'YData', 0);
                return;
            end

            minIc = min(min(Ic));
            maxIc = max(max(Ic));
            Icn = (Ic-minIc)/(maxIc-minIc);

            % Convert grid to black and white
            Nb = 10;
            [counts, centers] = hist(reshape(Icn,1,[]),Nb);
            % Filter some noise
            countsl = max(log(counts),0);
            countsf = (countsl(1:end-1)+countsl(2:end))/2;
            centersf = (centers(1:end-1)+centers(2:end))/2;
            
            % Find troughs between peaks
            med = median(countsf);
            vals = countsf;
            vals(vals>med) = 0;
            zones = [];
            nz1 = 0;
            while true
                nz0 = find(vals(nz1+1:end) ~= 0, 1, 'first');
                if isempty(nz0)
                    break
                else
                    nz0 = nz1 + nz0;
                end
                nz1 = find(vals(nz0+1:end) == 0, 1, 'first');
                if isempty(nz1)
                    break
                else
                    nz1 = nz0 + nz1;
                end
                zones = [zones;[nz0, nz1]]; %#ok<AGROW>
            end
            % Find largest trough
            len = diff(zones,[],2);
            largest = find(len == max(len), 1, 'first');
            if isempty(largest)
                I(grid_bb(1):grid_bb(3), grid_bb(2):grid_bb(4)) = Ic;
                set(WGV.ih,'CData',I);
                return;
            end
            thresh = 1.2*centersf(floor(mean(zones(largest,:))));
                        
            % Prepare plot data
            px = WGV.normalize(centersf,WGV.pxlim);
            py = WGV.normalize(countsf,WGV.pylim);
            set(WGV.ph, 'XData', px, 'YData', py);

            % Prepare preview image
            I = 0.25*grayImage;
            if thresh<0 || thresh>1
                I(grid_bb(1):grid_bb(3), grid_bb(2):grid_bb(4)) = Ic;
                set(WGV.ih,'CData',I);
                return;
            end

            Ibw = im2bw(Icn, thresh);

            % Remove tiles to leave only letters
            tiles = bwareaopen(~Ibw,floor(0.1*numel(Ibw)));
            Ilet = ~Ibw-tiles;
            % Clean up "debris"
            Ilet = bwareaopen(Ilet,80);
            
            I(grid_bb(1):grid_bb(3), grid_bb(2):grid_bb(4)) = 255*Ilet;

            % Store the image
            WGV.I_store{WGV.I_n} = Ilet;
            WGV.I_n = WGV.I_n + 1;
            if WGV.I_n > WGV.n_store
                WGV.I_n = 1;
            end

            % Compare images
            if ~any(cellfun(@isempty,WGV.I_store))
                c = zeros(1,WGV.n_store);
                for i = 1:WGV.n_store-1
                    c(i) = corr2(WGV.I_store{i}, ...
                        imresize(WGV.I_store{i+1},size(WGV.I_store{i})));
                end
                c(end) = corr2(WGV.I_store{end}, ...
                    imresize(WGV.I_store{1},size(WGV.I_store{end})));

                if all(c>0.8)
                    WGV.grid_stable = true;
                    
                    % Display image with a green border
                    h = floor(size(I,1)/20);
                    w = floor(size(I,2)/20);
                    I = double(repmat(I,[1,1,3]))/255;
                    
                    green = ones(size(I,1), size(I,2), 3);
                    green(:,:,[1, 3]) = 0.3;
                    green(h:end-h, w:end-w, :) = ...
                        I(h:end-h, w:end-w, :);
                    I = green;
                    
                    % Perform post processing to get individual letters
                    im_letters = WGV.post_process(Ilet);
                    % Get the words length
                    try
                        [im_word_len, im_word] = WGV.word_len_proc(Ic2);
                    catch ME
                        disp(['ERROR! ',ME.message])
                        disp(ME.cause)
                        for i=1:length(ME.stack)
                            disp(['In file: ',ME.stack(i).file,...
                                ', name: ',ME.stack(i).name,...
                                ', line: ',num2str(ME.stack(i).line)])
                        end
                        
                        im_word_len = [];
                        im_word = Ic2;
                    end
                    I(round(min_y):round(max_y), ...
                        grid_bb(2):grid_bb(4),:) = ...
                        repmat(double(im_word),[1,1,3]);
                    
                    if size(im_letters,2) == sum(im_word_len)
                        WGV.letters = im_letters;
                        WGV.word_len = im_word_len;
                        
                        % Store image for easy comparison with results
                        WGV.I_store{1} = [Ic;Ic2];

                        if ~isempty(WGV.capture_cb)
                            WGV.capture_cb(WGV);
                        end
                    else
                        WGV.grid_stable = false;
                    end
                else
                    WGV.grid_stable = false;
                end
            end

            % Display the image.
            set(WGV.ih,'CData',I);
            drawnow
        end
        
        function mserStats = GetStats(WGV, mserRegions, sz) %#ok<INUSL>
            % First, convert the x,y pixel location data within mserRegions into linear
            % indices as required by regionprops.
            pixelIdxList = cellfun(@(xy)sub2ind(sz, xy(:,2), xy(:,1)), ...
                mserRegions.PixelList, 'UniformOutput', false);

            % Next, pack the data into a connected component struct.
            mserConnComp.Connectivity = 8;
            mserConnComp.ImageSize = sz;
            mserConnComp.NumObjects = mserRegions.Count;
            mserConnComp.PixelIdxList = pixelIdxList;

            % Use regionprops to measure MSER properties
            mserStats = regionprops(mserConnComp, ...
                'Area', 'BoundingBox', 'Orientation', ...
                'Eccentricity', 'Solidity');
        end
        
        function im_letters = post_process(WGV, I)            
            % Get letters properties
            CC = regionprops(I,'Area','BoundingBox','PixelList');
            if length(CC)<4
                im_letters = [];
                return
            end
            
            % Find bounding box for each letter
            BB = reshape([CC.BoundingBox],4,[])';
            max_size = max(BB(:,[4,3]),[],1);

            % Get each letter into a standardized frame and resize it to
            % the desired size for NN
            proc_ltrs = zeros(WGV.des_size(1), WGV.des_size(2), length(CC));
            letter_pos = zeros(length(CC),2);
            for i = 1:length(CC)
                this_char = zeros(max_size(1), max_size(2));

                % find bounding box
                this_BB = [min(CC(i).PixelList),max(CC(i).PixelList)];
                Iregion = I(this_BB(2):this_BB(4), this_BB(1):this_BB(3));
                reg_size = size(Iregion);
                letter_pos(i,:) = [this_BB(2), this_BB(1)];

                % Position character
                delta = floor((max_size-reg_size)/2);
                this_char(1+delta(1):delta(1)+reg_size(1), ...
                    1+delta(2):delta(2)+reg_size(2)) = Iregion;

                % Scale to standard size (for NN)
                proc_ltrs(:,:,i) = ...
                    imresize(this_char, WGV.des_size, 'nearest');
            end
            
            % Position letters on grid
            r = ceil(sqrt(length(CC)));
            grid_min = min(letter_pos);
            grid_max = max(letter_pos);
            x_pos = linspace(grid_min(1), grid_max(1), r);
            y_pos = linspace(grid_min(2), grid_max(2), r);
            letter_sub = zeros(length(CC),2);
            for i = 1:length(CC)
                dist_y = abs(letter_pos(i,2)-y_pos);
                dist_x = abs(letter_pos(i,1)-x_pos);
                letter_sub(i,:) = ...
                    [find(dist_y == min(dist_y), 1, 'first'), ...
                     find(dist_x == min(dist_x), 1, 'first')];
            end

            % Prepare output vector
            im_letters = zeros(prod(WGV.des_size), length(CC));
            for i = 1:length(CC)
                ind = sub2ind([r,r],letter_sub(i,1),letter_sub(i,2));
                im_letters(:,ind) = reshape(proc_ltrs(:,:,i),[],1);
            end
        end
        
        function done_capture(WGV, ~, ~)
            if WGV.laptop
                closepreview(WGV.vid);
                delete(WGV.fh)
            else
                WGV.go_cap = false;
            end 
            
            if ~isempty(WGV.done_cb)
                WGV.done_cb(WGV)
            end
        end
    end
    
    methods(Static)
        function key_press_cb(~, ~, ~)
            k = get (gcf, 'CurrentKey');

            if strcmp(k,'return')
                % Save image to disk
                imwrite(rgbImage,'sample.png');
            end
        end
        
        function [word_len, Iout] = word_len_proc(I)
            minI = min(min(I));
            maxI = max(max(I));
            In = (I-minI)/(maxI-minI);
        %     figure; imshow(In);

            % Apply Canny edge detection
            [~, thresh] = edge(In,'canny');
            BW = edge(In,'canny',thresh*1.5);
        %     figure; imshow(BW)

            % Dilate
            se90 = strel('line', 5, 90);
            se0 = strel('line', 5, 0);
            BWdil = imdilate(BW, [se90 se0]);
            BWdil = ~bwareaopen(~BWdil,30);
        %     figure(); imshow(BWdil);

            % Clean "hint" tiles
            CC1 = bwconncomp(BWdil);
            if CC1.NumObjects > 1
                BB1 = zeros(CC1.NumObjects,4);
                for i = 1:CC1.NumObjects
                    [y,x] = ind2sub(CC1.ImageSize,CC1.PixelIdxList{i});
                    BB1(i,:) = [min(x), max(x), min(y), max(y)];
                end

                CC2 = bwconncomp(~BWdil);
                BB2 = zeros(CC2.NumObjects,4);
                for i = 1:CC2.NumObjects
                    [y,x] = ind2sub(CC2.ImageSize,CC2.PixelIdxList{i});
                    BB2(i,:) = [min(x), max(x), min(y), max(y)];
                end

                for i = 1:CC2.NumObjects
                    % Check if region is within a region in ~BWdil
                    inside = zeros(1,CC1.NumObjects);
                    for j = 1:CC1.NumObjects
                        if BB2(i,1)<=BB1(j,1) && BB1(j,2)<=BB2(i,2) ...
                                && BB2(i,3)<=BB1(j,3) && BB1(j,4)<=BB2(i,4)
                            inside(j) = 1;
                        end
                    end
                    IDs = find(inside == 1);
                    if length(IDs) == 1
                        % Found it!
                        BWdil(CC2.PixelIdxList{i}) = 1;
                        BWdil(CC1.PixelIdxList{IDs}) = 0;
            %             figure; imshow(BWdil)
            %             figure; imshow(~BWdil(BB2(i,3):BB2(i,4), ...
            %                                   BB2(i,1):BB2(i,2)))
            %             figure; imshow(BWdil(BB1(IDs,3):BB1(IDs,4), ...
            %                                  BB1(IDs,1):BB1(IDs,2)))
                    end
                end
            end
    
            % Get connected components
            CC = bwconncomp(~BWdil);
            area = cellfun(@length,CC.PixelIdxList);
            Iout = ~BWdil;
            Iout(CC.PixelIdxList{area == max(area)}) = 0;
            CC.PixelIdxList(area == max(area)) = [];
            CC.NumObjects = length(CC.PixelIdxList);

            % Remove outliers
            prop = regionprops(CC,'Solidity','Eccentricity');
            sol = [prop.Solidity];
            ecc = [prop.Eccentricity];
            CC.PixelIdxList(sol < 0.75 | ecc > 0.8) = [];
            CC.NumObjects = length(CC.PixelIdxList);

            % Remove bridges (space between words)
            area = cellfun(@length,CC.PixelIdxList);
            area_mean = mean(area);
            area_std = std(area);
            CC.PixelIdxList(area > area_mean+3*area_std | ...
                area < area_mean-3*area_std) = [];

            % Reconstruct image with boxes only
            Ic = zeros(CC.ImageSize);
            for i = 1:length(CC.PixelIdxList)
                Ic(CC.PixelIdxList{i}) = 1;
            end
%             figure; imshow(Ic)
            
            % Join boxes horizontally to create words
            se0 = strel('line', 10, 0);
            Icd = imdilate(Ic,[se0 se0]);

            % Get words length
            CC = bwconncomp(Icd);
            props = regionprops(CC,'BoundingBox');
            N_words = length(props);
            bb = vertcat(props.BoundingBox);
            bb = [floor(bb(:,1:2)) floor(bb(:,3:4))];
            bb = max(bb,1);

            if N_words == 1
                word_ind = 1;
            else
                word_pos = [bb(:,1:2), (1:N_words)'];
                miny = min(word_pos(:,2));
                maxy = max(word_pos(:,2));
                distx = max(word_pos(:,1))-min(word_pos(:,1));
                r = ceil(sqrt(N_words));
                max_dist = (maxy-miny)/r/2;
                word_ind = zeros(1, N_words);

                if max_dist/distx < 0.01
                    % Single row
                    row = sortrows(word_pos,1);
                    word_ind = row(:,end);
                else
                    w = 1;
                    while ~isempty(word_pos)
                        % Get all words in row
                        dist = word_pos(:,2) - miny;
                        ids = find(abs(dist)<max_dist);
                        if isempty(ids)
                            warning('word_len_proc: No words found on row')
                            break
                        end
                        row = sortrows(word_pos(ids, [1 3]),1);
                        word_ind(w:w+size(row,1)-1) = row(:,end);

                        w = w+size(row,1);
                        word_pos(ids, :) = [];
                        miny = min(word_pos(:,2));
                    end
                end
            end
            
            % Rearrange the regions
            bb = bb(word_ind,:);

            word_len = zeros(1, N_words);
            for i = 1:N_words
                try
                    Iw = Ic(bb(i,2):bb(i,2 )+bb(i,4), bb(i,1):bb(i,1)+bb(i,3));
                catch
                    disp(bb(i,:))
                    disp(['min: ', int2str(bb(i,2)), ' - max: ', ...
                        int2str(bb(i,2)+bb(i,4)), 'min: ', int2str(bb(i,1)), ...
                        ' - max: ',int2str(bb(i,1)+bb(i,3))])
                    disp(size(Ic))
                end
                CC = bwconncomp(Iw);
                word_len(i) = length(CC.PixelIdxList);
            end
            if any(word_len<2)
                word_len = [];
            end
    
%             % Apply fuzzy logic edge detection
%             Gx = [-1 1];
%             Gy = Gx';
%             Ix = conv2(In,Gx,'same');
%             Iy = conv2(In,Gy,'same');
%             edgeFIS = newfis('edgeDetection');
%             edgeFIS = addvar(edgeFIS,'input','Ix',[-1 1]);
%             edgeFIS = addvar(edgeFIS,'input','Iy',[-1 1]);
% 
%             sx = 0.05; sy = 0.05;
%             % sx and sy specify the standard deviation for the zero
%             % membership function for the Ix and Iy inputs. You can
%             % change the values of sx and sy to adjust the edge
%             % detector performance. Increasing the values makes the
%             % algorithm less sensitive to the edges in the image and
%             % decreases the intensity of the detected edges.
%             edgeFIS = addmf(edgeFIS,'input',1,'zero','gaussmf',[sx 0]);
%             edgeFIS = addmf(edgeFIS,'input',2,'zero','gaussmf',[sy 0]);
%             edgeFIS = addvar(edgeFIS,'output','Iout',[0 1]);
% 
%             wa = 0.1; wb = 1; wc = 1;
%             ba = 0; bb = 0; bc = .7;
%             % The triplets specify the start, peak, and end of the
%             % triangles of the membership functions.
%             edgeFIS = addmf(edgeFIS,'output',1,'white','trimf',[wa wb wc]);
%             edgeFIS = addmf(edgeFIS,'output',1,'black','trimf',[ba bb bc]);
% 
%             r1 = 'If Ix is zero and Iy is zero then Iout is white';
%             r2 = 'If Ix is not zero or Iy is not zero then Iout is black';
%             r = char(r1,r2);
%             edgeFIS = parsrule(edgeFIS,r);
% 
%             Ieval = zeros(size(In));% Preallocate the output matrix
%             for ii = 1:size(In,1)
%                 Ieval(ii,:) = evalfis([(Ix(ii,:));(Iy(ii,:));]',edgeFIS);
%             end
% 
%             % figure()
%             % imshow(Ieval)
% 
%             % Transform to black and white
%             BWs = ~im2bw(Ieval, mean(reshape(Ieval,1,[])));
%             BWs = bwareaopen(BWs,10);
%             % figure()
%             % imshow(BWs)
% 
%             % Dilate
%             se90 = strel('line', 4, 90);
%             se0 = strel('line', 4, 0);
%             BWsdil = imdilate(BWs, [se90 se0]);
%             BWsdil = ~bwareaopen(~BWsdil,10);
%             % figure()
%             % imshow(BWsdil)
% 
% %             % Erode
% %             se90 = strel('line', 6, 90);
% %             se0 = strel('line', 6, 0);
% %             BWsero = imerode(BWsdil, [se90 se0]);
% %             % figure()
% %             % imshow(BWsero)
% 
%             % Get connected components
%             CC = bwconncomp(~BWsdil);
%             area = cellfun(@length,CC.PixelIdxList);
%             Iout = ~BWsdil;
%             Iout(CC.PixelIdxList{area == max(area)}) = 0;
%             CC.PixelIdxList(area == max(area)) = [];
% 
%             % Remove bridges
%             area = cellfun(@length,CC.PixelIdxList);
%             area_mean = mean(area);
%             area_std = std(area);
%             CC.PixelIdxList(area > area_mean+3*area_std) = [];
% 
%             % Reconstruct image with boxes only
%             Ic = zeros(CC.ImageSize);
%             for i = 1:length(CC.PixelIdxList)
%                 Ic(CC.PixelIdxList{i}) = 1;
%             end
        end
        
        function ndata = normalize(data, lim)
            % Transform to 0-1 range
            ndata = (data-min(data))/(max(data)-min(data));
            % Transform to lim range
            ndata = lim(1)+(lim(2)-lim(1))*ndata;
        end
        
        function thresh = imthresh(I)
            Nb = 10;
            [counts, centers] = hist(reshape(I,1,[]),Nb);
            % Filter some noise
            countsl = max(log(counts),0);
            countsf = (countsl(1:end-1)+countsl(2:end))/2;
            centersf = (centers(1:end-1)+centers(2:end))/2;
            
            % Find troughs between peaks
            med = median(countsf);
            vals = countsf;
            vals(vals>med) = 0;
            zones = [];
            nz1 = 0;
            while true
                nz0 = find(vals(nz1+1:end) ~= 0, 1, 'first');
                if isempty(nz0)
                    break
                else
                    nz0 = nz1 + nz0;
                end
                nz1 = find(vals(nz0+1:end) == 0, 1, 'first');
                if isempty(nz1)
                    break
                else
                    nz1 = nz0 + nz1;
                end
                zones = [zones;[nz0, nz1]]; %#ok<AGROW>
            end
            % Find largest trough
            len = diff(zones,[],2);
            largest = find(len == max(len), 1, 'first');
            if isempty(largest)
                thresh = [];
                return;
            end
            thresh = 1.2*centersf(floor(mean(zones(largest,:))));
        end
    end
end

