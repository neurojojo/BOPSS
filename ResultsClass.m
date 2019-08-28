classdef ResultsClass < handle
    
       properties
          ID
          Images
          Name
          Cells
          Puncta
          DetectedPuncta
          NegativeControl=[]
          ResultsTable=[]
          % Parameters for the analysis
          minsize
          maxsize
          celldilation
          regularResult
          decodedResult
          trainFrac=0.5
       end
       
       properties (Access=private)
           Pct_Cutoff=95
           B
           Bdecoded
           L
           Nimages
           include
           analyzerData
           grayImages
           showDistributions=0
           CellAreas_forML
           PunctaAreas_forML
           resultsTableForExport
           FilesTable
       end
       
       methods
           
           function obj = ResultsClass(varargin) % Constructor Method
               folder = varargin{1};
               files = dir(folder);
               counter=1;
               datatype = {'Data','Negative Control'};
    
               for i = 1:size(files,1)
                   if ~isempty(strfind(files(i).name,'.tif')) | ~isempty(strfind(files(i).name,'.Tiff'))
                       this_image = imread([folder,'/',files(i).name]);
                       obj.Images{counter,1} = this_image(:,:,[1:3]);
                       obj.Name{counter} = files(i).name; % MATLAB2016a needs a cell here for Name assignment
                       if ~isempty(strfind(files(i).name,'neg')) | ~isempty(strfind(files(i).name,'Neg'))
                           obj.NegativeControl(counter) = 2; else; obj.NegativeControl(counter) = 1;
                       end
                       counter=counter+1;
                   end
               end
               
               %result.FilesTable = table([1:size(results.Images,1)]',results.Name',datatype(results.NegativeControl)');
               obj.ID = sprintf('results_%s',datestr(datetime,'ddmmm_HHMM'));
           end  % End of constructor Method
           
           % This method accomplishes the color-based segmentation
           function output = std_(obj,input)
               output = std(double(input));
           end
           
           function obj=thresholdImages(obj)
              cform = makecform('srgb2lab');
                for i = 1:size(obj.Images,1)
                    new_cform = []; % This is important because it allows images to be different sizes
                    im1_cform = double(applycform(obj.Images{i},cform));
                    im1_cform_comp = double(imcomplement(im1_cform(:,:,1)));
                    im1_thresh = mat2gray(im1_cform_comp,[min(im1_cform_comp(:)),max(im1_cform_comp(:))]);
                    
                    threshold = graythresh(im1_thresh);
                    im1_thresh = imfill(double(im2bw(im1_thresh,threshold)),'holes');
                    im1_thresh(im1_thresh==0)=NaN;
                    
                    new_cform(:,:,1) = im1_thresh.*im1_cform(:,:,2);
                    new_cform(:,:,2) = im1_thresh.*im1_cform(:,:,3);
                    
                    q1 = obj.flatmat(new_cform(:,:,1));
                    q2 = obj.flatmat(new_cform(:,:,2));
                    
                    %q3 = obj.flatmat(im1_thresh.*imfilter(im1_cform(:,:,2),fspecial('gaussian',20,10)));
                    %q4 = obj.flatmat(im1_thresh.*imfilter(im1_cform(:,:,3),fspecial('gaussian',20,10)));
                    
                    warning('off','stats:kmeans:MissingDataRemoved');
                    [sorted,idx] = kmeans([q1,q2],2,'replicates',10);
                    warning('on','stats:kmeans:MissingDataRemoved');
                    
                    tmp1 = reshape(sorted,size(obj.Images{i},1),size(obj.Images{i},2)); 
                    tmp2 = reshape(sorted,size(obj.Images{i},1),size(obj.Images{i},2)); 

                    tmp1(tmp1~=1)=0; % Two groups for K-means, set one equal to one
                    tmp2(tmp2~=2)=0; % Other equal to two, we don't yet know the ID of groups
                    tmp2(tmp2==2)=1;

                    tmp = zeros(size(obj.Images{i},1),size(obj.Images{i},2));
                    tmp(:,:,1) = tmp1;
                    tmp(:,:,2) = tmp2;

                    % Discern puncta from cells by mean area
                    [~,tmp1_L] = bwboundaries(tmp1,'noholes');
                    [~,tmp2_L] = bwboundaries(tmp2,'noholes');

                    tmp_area{1} = cell2mat(struct2cell(regionprops(tmp1_L,'Area')));
                    tmp_area{2} = cell2mat(struct2cell(regionprops(tmp2_L,'Area')));
                    
                    if (idx(1,1)<idx(2,1)) & (idx(1,2)>idx(2,2))
                        puncta = 1; cells = 2;
                    else
                        puncta = 2; cells = 1;
                    end
    
                    obj.Cells{i} = tmp(:,:,cells);
                    obj.Puncta{i} = tmp(:,:,puncta);
                    
                    obj.CellAreas_forML = [obj.CellAreas_forML,tmp_area{cells}];
                    obj.PunctaAreas_forML = [obj.PunctaAreas_forML,tmp_area{puncta}];
                    
                    % Create grayscale images for later use with
                    % optimization if desired
                    obj.grayImages{i} = rgb2gray(obj.Images{i});
                end
           end
           
           
           % This method uses the input parameters to accomplish the
           % thresholding
           function obj=countPuncta(obj)
               
                for i = 1:size(obj.Images,1)
                    
                    Cells = obj.Cells{i};
                    Puncta = obj.Puncta{i};
                    
                    % Automatically infer cell dilation and puncta max size based on the
                    % distribution from the cell and puncta images
                    obj.celldilation = sqrt(mean(obj.CellAreas_forML(obj.CellAreas_forML>prctile(obj.PunctaAreas_forML,95))));
                    obj.maxsize = fix(prctile(obj.PunctaAreas_forML,obj.Pct_Cutoff));
                    
                    dilate_each_cell = fix(obj.celldilation);
                    
                    dilated_cells = bwareaopen(Cells,obj.maxsize)+bwareaopen(Puncta,obj.maxsize);
                    dilated_cells = imfill(imdilate(dilated_cells,strel('disk',dilate_each_cell)),'holes');
                    
                    if isempty(obj.minsize); obj.minsize=2; % Cannot allow single pixel puncta
                    end
                    
                    cell_less_image = Puncta.*imcomplement(dilated_cells);

                    [B,L] = bwboundaries(cell_less_image,'noholes');
                    puncta_area = cell2mat(struct2cell(regionprops(L,'area')));

                    include = setxor([1:numel(B)],[find(puncta_area>obj.maxsize),...
                        find(puncta_area<obj.minsize)]);
                    obj.DetectedPuncta(i,1) = numel(include);
                    
                    % Private properties to be used for optimization if
                    % desired
                    obj.B{i} = B(include);
                    obj.L{i} = ismember(L,include); % Only extracts selected puncta
                    obj.include{i} = include;
                    
                end
                % Line 118 is fix for MATLAB 2016a
                if size(obj.Name)==fliplr(size(obj.DetectedPuncta)); obj.DetectedPuncta = obj.DetectedPuncta'; end
                T = table(obj.Name,obj.DetectedPuncta);
                T.Properties.VariableNames = {'Filename',sprintf('Detected_puncta')};
                obj.ResultsTable = T;
                fprintf('\nAnalyzed using inferred puncta size: %1.0f, max puncta size: %1.0f, cell dilation: %1.0f\n',obj.minsize,obj.maxsize,obj.celldilation);
               
           end
            
           function analyze(obj,varargin)
               
               % If there is no specification of which files are positive
               % control and which are negative control, this _if_ 
               % will fill in indices automatically
               if nargin==1; fp_indices = find(obj.NegativeControl==2); tp_indices=find(obj.NegativeControl==1); end
               if isempty(obj.Cells); obj = obj.thresholdImages(); obj = obj.countPuncta(); end
               if isempty(obj.ResultsTable); obj = obj.countPuncta(); end
               
               % Create two matrices containing statistics about the puncta
               analyzerData_color = []; % This matrix contains statistics about the color values of puncta
               analyzerData_bw = []; % This matrix contains statistics about the shape of the puncta irrespective of intensity 
               obj.Nimages = size(obj.grayImages,2);
               
               % Getting information about identified objects
               for i = 1:obj.Nimages
                   shape_tmp_struct = ((regionprops(obj.L{i},...
                            'convexarea',...        % Property 5
                            'eccentricity',...      % Property 6
                            'majoraxislength',...   % Property 7
                            'minoraxislength',...   % Property 8
                            'orientation')));   
                   shape_tmp = [];
                   shape_tmp(1,:) = vertcat( shape_tmp_struct.ConvexArea );
                   shape_tmp(2,:) = vertcat( shape_tmp_struct.Eccentricity );
                   shape_tmp(3,:) = vertcat( shape_tmp_struct.MajorAxisLength );
                   shape_tmp(4,:) = vertcat( shape_tmp_struct.MinorAxisLength );
                   shape_tmp(5,:) = vertcat( shape_tmp_struct.Orientation );
                   
                   loc_tmp = regionprops(obj.L{i},'Centroid'); % Locations of all puncta
                   loc_tmp = reshape(cell2mat(struct2cell(loc_tmp)),2,size(shape_tmp,2));
                   
                   % Compute pairwise distances
                   pdistmatrix = squareform(pdist([loc_tmp(1,:)',loc_tmp(2,:)']));
                   c = size(pdistmatrix, 1); idx = 1:c+1:numel(pdistmatrix); pdistmatrix(idx)=NaN;
                   
                   % Compute # of partners within maxsize of a puncta
                   tmp = [shape_tmp; ... % Properties 5-9
                       log(min(pdistmatrix,[],1)); ... % Property 10
                       sum(pdistmatrix<obj.maxsize)]; % Property 11
                   analyzerData_bw = [analyzerData_bw, tmp];
                   
                   pixel_info{i} = regionprops(obj.L{i},obj.grayImages{i},'pixelvalues');
                   
               end
               
               
               for i = 1:obj.Nimages
                   % We parse out three channels for the data
                   tmp=struct2cell(pixel_info{i}); % Looks at the pixel values for each image
                   
                   [punctamean,punctasd,punctaentropy,punctarange] = deal(cellfun(@mean,tmp),cellfun(@obj.std_,tmp),cellfun(@entropy,tmp),cellfun(@range,tmp));
                   analyzerData_color = cat(2,analyzerData_color,...
                       [repmat(i,1,numel(punctamean));... % Labels for the image index in the sequence
                       double(punctamean);...
                       (double(punctasd));...
                       double(punctaentropy);...
                       (double(punctarange))]);
               end
               
               analyzerData_bw(1,:) = (analyzerData_bw(1,:));
               analyzerData_bw(2,:) = (analyzerData_bw(2,:));
               analyzerData_bw(5,:) = (analyzerData_bw(5,:));
               analyzerData = [analyzerData_color;analyzerData_bw];
               Nrois = size(analyzerData,2);
               Nfeatures = size(analyzerData,1);
               
               % LEARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               fp = find( ismember(analyzerData(1,:),fp_indices) == 1); 
               tp = find( ismember(analyzerData(1,:),tp_indices) == 1); 
               
               % Subsample the false positives and true positives
               fp_train = fp( randi( numel(fp), fix(obj.trainFrac*numel(fp)), 1) );
               tp_train = tp ( randi( numel(tp), fix(obj.trainFrac*numel(tp)), 1) );
               
               tmp = analyzerData([2:end],[fp_train,tp_train]);
               tmp_rs = (tmp-min(tmp,[],2))./range(tmp,2);
               
               % Extract these 
               fp_analyzer = analyzerData(:,fp_train);
               tp_analyzer = analyzerData(:,tp_train);
               
               fp_means = mean(fp_analyzer,2);
               tp_means = mean(tp_analyzer,2);
               
               fp_sds = std(double(fp_analyzer),[],2);
               tp_sds = std(double(tp_analyzer),[],2);
               
               for i = 1:Nfeatures % Attributes (starting at row 3)
                   attribute_bins{i} = [min(min(analyzerData(i,:))) - .1*iqr(analyzerData(i,:))*numel(analyzerData(i,:))^(-1/3):2*iqr(analyzerData(i,:))*numel(analyzerData(i,:))^(-1/3):1.2*max(max(analyzerData(i,:)))];
               end
               
               % END LEARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           % OPTIONAL PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if obj.showDistributions; figure;plt_titles = {'id','punctamean','log(punctasd)','punctaentropy','punctarange','convexarea','eccentricity','majoraxislength','minoraxislength','orientation','nnd','numnears'};
               for i = 1:size(analyzerData,1);
                       H_tp = (1/sum(tp))*histc(analyzerData(i,tp_train),attribute_bins{i});
                       H_fp = (1/sum(fp))*histc(analyzerData(i,fp_train),attribute_bins{i});
                       tmp = corrcoef(H_tp,H_fp); cc(i)=tmp(2);
                       subplot(3,4,i); plot(attribute_bins{i},H_tp,'-'); hold on; plot(attribute_bins{i},H_fp,'-'); hold on;
                       title(sprintf('%1.2f %s', cc(i), plt_titles{i}));
                       box off;set(gca,'YTick',[],'YColor','w');
               end
           end
           % END OPTIONAL PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           % Run the analysis on *all* ROIs (train and test -- for
           % visualization purposes)
           %
           % However, for quantification, only test data is placed into output

           for i = 1:Nfeatures % Attributes
               for j = 1:Nrois % Observations
                   output1(i,j) = normpdf(analyzerData(i,j),tp_means(i),tp_sds(i))./ ...
                           sum(normpdf(attribute_bins{i},tp_means(i),tp_sds(i)));
                   output2(i,j) = normpdf(analyzerData(i,j),fp_means(i),fp_sds(i))./ ...
                           sum(normpdf(attribute_bins{i},fp_means(i),fp_sds(i)));
               end
           end

            output1 = sum(log(output1),1);
            output2 = sum(log(output2),1);
            decoded = (output1>output2).*analyzerData(1,:);
            
            test_idxes = setxor([1:Nrois],[tp_train,fp_train]);
            decoded_on_test = decoded(test_idxes);
            
            for i = 1:obj.Nimages; obj.Bdecoded{i} = obj.B{i}(gt(decoded(1,analyzerData(1,:)==i),0)); end
                
            obj.regularResult = histc(analyzerData(1,test_idxes),[1:obj.Nimages]);
            obj.decodedResult = histc(decoded_on_test(decoded_on_test>0),[1:obj.Nimages]);
            
          obj.resultsTableForExport = table(obj.Name',obj.regularResult',obj.decodedResult','variablenames',{'Name','Undecoded','Decoded'});
          disp(obj.resultsTableForExport);
                   
                   
           end
           
           function saveImages(obj,varargin)
                    
               if sum(strcmp(varargin,'decoded')); B = obj.Bdecoded; analysistype = 'decoded'; else; B = obj.B; analysistype = 'undecoded'; end
               if sum(strcmp(varargin,'show')); show = 'on'; else; show = 'off'; end
                    for i = 1:size(obj.Name,2)
                        thisB = B{i};
                        figure('visible',show); imagesc(obj.Images{i}); 
                        hold on; 
                        for j = 1:size(thisB,1)
                            plot(thisB{j}(:,2),thisB{j}(:,1),'r-','linewidth',.5); 
                        end
                        title(strrep(sprintf('%s_%s',analysistype,obj.Name{i}),'_',''));
                        print(sprintf('%s_%s',analysistype,obj.Name{i}),'-dtiff');
                        fprintf('File saved to %s_%s\n',analysistype,obj.Name{i});
                    end
           end
           
           function makeFigure(obj)
                figure('color','w');
                set(gcf,'position',[209,147,1156,615]);
                subplot(1,2,1); bar(obj.regularResult); title('Undecoded result')
                ylabel('Number of puncta detected')
                set(gca,'XTickLabel',strrep(strrep(obj.Name,'_',''),'.tif','')); 
                set(gca,'XTickLabelRotation',45)
                subplot(1,2,2); bar(obj.decodedResult); title('Decoded result')
                ylabel('Number of puncta detected')
                set(gca,'XTickLabel',strrep(strrep(obj.Name,'_',''),'.tif','')); 
                set(gca,'XTickLabelRotation',45)
           end
           
           function saveResults(obj,filename)
               if nargin==1; filename = obj.ID; end
               writetable(obj.resultsTableForExport,filename)
           end
           
           % Auxillary methods for accomplishing other things
           function out = flatmat(obj,in)
              out = in(:); 
           end
           
       end % Of all methods
end