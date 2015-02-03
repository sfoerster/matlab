opts = [];
opts.importData = true;
opts.graphTemp = false;
opts.graphTempN = false;
opts.trainSVM = true;
opts.classTree = true;
opts.boostTree = false;
opts.trainANN = true;
opts.plotANN = true;

params = [];

%SVM
params.tdset = 'TAVG';
params.tempN = 1759; %3026;
params.latRange = [33 55]; %[33 64];
params.lonRange = [-131 -100]; %[-131 -70];
params.tempThresh = 4;
params.svmKernel = 'linear'; %'rbf';

%ANN
params.ANNnodes = 20;


% files = {{'TAVG','Complete_TAVG_EqualArea.nc'}, ...
%     {'TMAX','Complete_TMAX_EqualArea.nc'}, ...
%     {'TMIN','Complete_TMIN_EqualArea.nc'}};
files = {{'TAVG','Complete_TAVG_EqualArea.nc'}};


if opts.importData
    DATA = struct([]);
    
    % import temperature data
    for n=1:length(files)

        t_info = ncinfo(files{n}{2});

        for m=1:length(t_info.Variables)
            DATA(1).(files{n}{1})(1).(t_info.Variables(1,m).Name) = ...
                double(ncread(files{n}{2},t_info.Variables(1,m).Name));

            if strcmp(t_info.Variables(1,m).Name,'time')
                tmp = DATA(1).(files{n}{1})(1).(t_info.Variables(1,m).Name);
                DATA(1).(files{n}{1})(1).(t_info.Variables(1,m).Name) = ...
                    doy2date((tmp-floor(tmp))*365,floor(tmp));
            end
        end

    end
    
    % import sonar data
    DATA(1).sonar(1) = load('supervised_data.mat');
    
    % import ionosphere data
    DATA(1).ionosphere(1) = load('ionosphere.mat');
end
    
if opts.graphTemp
    for n=1:size(DATA(1).(params.tdset)(1).temperature,2)
        figure(23); clf;
        scatter3(DATA(1).(params.tdset)(1).longitude, ...
            DATA(1).(params.tdset)(1).latitude, ...
            DATA(1).(params.tdset)(1).temperature(:,n), ...
            5, ...
            DATA(1).(params.tdset)(1).temperature(:,n)); 

        view(0,90); 
        ylabel('latitude'); 
        xlabel('longitude'); 
        tmp = DATA(1).(params.tdset)(1).time(n,1);
        title({['Time: ' datestr(tmp)],['n: ' num2str(n)]})

    end
end

if opts.graphTempN
    
    n = params.tempN;
    
    figure(23); clf;
    scatter3(DATA(1).(params.tdset)(1).longitude, ...
        DATA(1).(params.tdset)(1).latitude, ...
        DATA(1).(params.tdset)(1).temperature(:,n), ...
        5, ...
        DATA(1).(params.tdset)(1).temperature(:,n), ...
        'o'); 

    view(0,90); 
    ylabel('\circ latitude'); 
    xlabel('\circ longitude'); 
    tmp = DATA(1).(params.tdset)(1).time(n,1);
    title({['Time: ' datestr(tmp)],['n: ' num2str(n)]})
    hcb = colorbar;
    cbtitle = get(hcb,'Title');
    set(cbtitle,'String','\circ Celsius')

end

if opts.trainSVM
    
    n = params.tempN;
    
    lon = DATA(1).(params.tdset)(1).longitude;
    lat = DATA(1).(params.tdset)(1).latitude;
    temp = DATA(1).(params.tdset)(1).temperature(:,n);
    
    lonIdx = lon > params.lonRange(1) & lon < params.lonRange(2);
    latIdx = lat > params.latRange(1) & lat < params.latRange(2);
    
    lon = lon(lonIdx & latIdx);
    lat = lat(lonIdx & latIdx);
    temp = temp(lonIdx & latIdx);
    
%     scatter3(lon,lat,temp,5,temp);
%     view(0,90)
%     colormap(jet)
%     colorbar

    htemp = double(temp > params.tempThresh); % threshold at 5 degrees C
    htemp(htemp==0) = -1;
    data3 = [lon lat];
    
    %Train the SVM Classifier
    cl = fitcsvm(data3,htemp,'Standardize',true,'KernelFunction',params.svmKernel,...
        'KernelScale','auto','BoxConstraint',Inf,'ClassNames',[-1,1]);

    % Predict scores over the grid
    [x1Grid,x2Grid] = meshgrid(linspace(min(data3(:,1)),max(data3(:,1)),100),...
        linspace(min(data3(:,2)),max(data3(:,2)),100));
    xGrid = [x1Grid(:),x2Grid(:)];
    [~,scores] = predict(cl,xGrid);

    % Plot the data and the decision boundary
    figure;
    h(1:2) = gscatter(data3(:,1),data3(:,2),htemp,'rb','.');
    hold on
    ezpolar(@(x)1);
    h(3) = plot(data3(cl.IsSupportVector,1),data3(cl.IsSupportVector,2),'ko');
    contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
    legend(h,{'-1','+1','Support Vectors'});
    xlabel('\circ longitude')
    ylabel('\circ latitude')
    title({['Support Vector Machine around ' num2str(params.tempThresh) '\circ Celsius region'], ...
        ['Kernel: ' params.svmKernel], ...
        ['Date: ' datestr(DATA(1).(params.tdset)(1).time(n,1))]})
    axis equal
    hold off

    
end

if opts.classTree
    
   %load ionosphere % contains X and Y variables
   ctree1 = ClassificationTree.fit(DATA(1).ionosphere(1).X, ...
       DATA(1).ionosphere(1).Y); 
   view(ctree1,'mode','graph') 
   
   ctree2 = ClassificationTree.fit(DATA(1).ionosphere(1).X, ...
       DATA(1).ionosphere(1).Y, ...
       'minleaf',40);
   view(ctree2,'mode','graph')
   
   % find optimal pruning level
   [~,~,~,bestlevel] = cvLoss(ctree1, ...
    'SubTrees','All','TreeSize','min');

    ctree1p = prune(ctree1,'Level',bestlevel);
    view(ctree1p,'mode','graph')
    
    % GentleBoost
    gb = fitensemble(DATA(1).ionosphere(1).X, ...
       DATA(1).ionosphere(1).Y, ...
       'GentleBoost', ...
       300, ...
       'Tree');
    
   figure;
   plot(resubLoss(gb,'mode','cumulative'),'b-');
   xlabel('Number of trees');
   ylabel('Resubstitution error');
   grid on
   legend('GentleBoost','Location','NE');
   
   figure;
   gbcv = crossval(gb,'kfold',5);
   plot(kfoldLoss(gbcv,'mode','cumulative'),'b-');
   xlabel('Number of trees');
   ylabel('Cross-validated error');
   grid on
   legend('GentleBoost','Location','NE');
   
   sbspc = fitensemble(DATA(1).ionosphere(1).X, ...
       DATA(1).ionosphere(1).Y, ...
       'Subspace', ... 
       300, ...
       'KNN');
   
   figure;
   plot(resubLoss(sbspc,'mode','cumulative'),'b-');
   xlabel('K');
   ylabel('Resubstitution error');
   grid on
   legend('Subspace - KNN','Location','NE');
   
   figure;
   sbcv = crossval(sbspc,'kfold',5);
   plot(kfoldLoss(sbcv,'mode','cumulative'),'b-');
   xlabel('K');
   ylabel('Cross-validated error');
   grid on
   legend('Subspace - KNN','Location','NE');
    
end

if opts.boostTree
    
   
%     s = urlread(['http://archive.ics.uci.edu/ml/' ...
%        'machine-learning-databases/adult/adult.data']);
% 
%     %'adult.data' represents missing data as '?'. Replace instances of missing data with the blank string '':
%     s = strrep(s,'?','');
% 
%     %Put the data into a MATLAB dataset array:
%     fid = fopen('adult.txt','w');
%     fwrite(fid,s);
%     fclose(fid);
%     clear s;
    VarNames = {'age' 'workclass' 'fnlwgt' 'education' 'education_num' ...
        'marital_status' 'occupation' 'relationship' 'race' ...
        'sex' 'capital_gain' 'capital_loss' ...
        'hours_per_week' 'native_country' 'income'};
    ds = dataset('file','adult.txt','VarNames',VarNames,...
        'Delimiter',',','ReadVarNames',false,'Format',...
       '%u%s%u%s%u%s%s%s%s%s%u%u%u%s%s');
    cat = ~datasetfun(@isnumeric,ds(:,1:end-1)); % Logical indices
    %                                             of categorical variables
    catcol = find(cat); % indices of categorical variables

    %Many predictors in the data are categorical. Convert those fields in the dataset array to nominal:
    ds.workclass = nominal(ds.workclass);
    ds.education = nominal(ds.education);
    ds.marital_status = nominal(ds.marital_status);
    ds.occupation = nominal(ds.occupation);
    ds.relationship = nominal(ds.relationship);
    ds.race = nominal(ds.race);
    ds.sex = nominal(ds.sex);
    ds.native_country = nominal(ds.native_country);
    ds.income = nominal(ds.income);

    %Convert the dataset array into numerical variables for fitensemble:
    X = double(ds(:,1:end-1));
    Y = ds.income;

    %Some variables have many levels. Plot the number of levels of each predictor:
    ncat = zeros(1,numel(catcol));
    for c=1:numel(catcol)
        [~,gn] = grp2idx(X(:,catcol(c)));
        ncat(c) = numel(gn);
    end
    figure;
    bar(catcol,ncat);
    xlabel('Predictor');
    ylabel('Number of categories');

    %Predictor 14 ('native_country') has more than 40 categorical levels. 
    %This is too many levels for any method except LogitBoost and GentleBoost.

    %Create classification ensembles using both LogitBoost and GentleBoost:
    lb = fitensemble(X,Y,'LogitBoost',300,'Tree','CategoricalPredictors',cat,...
        'PredictorNames',VarNames(1:end-1),'ResponseName','income');
    gb = fitensemble(X,Y,'GentleBoost',300,'Tree','CategoricalPredictors',cat,...
        'PredictorNames',VarNames(1:end-1),'ResponseName','income');

    %Examine the resubstitution error for the two ensembles:
    figure;
    plot(resubLoss(lb,'mode','cumulative'));
    hold on
    plot(resubLoss(gb,'mode','cumulative'),'r--');
    hold off
    xlabel('Number of trees');
    ylabel('Resubstitution error');
    legend('LogitBoost','GentleBoost','Location','NE');

    %The algorithms have similar resubstitution error.

    %Estimate the generalization error for the two algorithms by cross validation.
    lbcv = crossval(lb,'kfold',5);
    gbcv = crossval(gb,'kfold',5);
    figure;
    plot(kfoldLoss(lbcv,'mode','cumulative'));
    hold on
    plot(kfoldLoss(gbcv,'mode','cumulative'),'r--');
    hold off
    xlabel('Number of trees');
    ylabel('Cross-validated error');
    legend('LogitBoost','GentleBoost','Location','NE'); 
    
end

if opts.trainANN
    
    % sonar
    net = feedforwardnet(params.ANNnodes,'trainlm');
    sonarInput = [DATA(1).sonar(1).NCC DATA(1).sonar(1).SINR];
    sonarTarget = DATA(1).sonar(1).LLHR;
    net = train(net, ...
        sonarInput', ...
        sonarTarget');
    view(net)
    y = net(sonarInput');
    perf = perform(net,y,sonarTarget');
      
    
end

if opts.plotANN
   
    [NCCfield,SINRfield] = meshgrid(linspace(min(sonarInput(:,1)),max(sonarInput(:,1)),100), ...
        linspace(min(sonarInput(:,2)),max(sonarInput(:,2)),100));
    
    tfield = reshape(net([NCCfield(:) SINRfield(:)]'),100,100);
    
    figure(27); clf;
    surf(NCCfield,SINRfield,tfield)
    xlabel('Normalized Cross-Correlation')
    ylabel('SINR')
    zlabel('LLHR')
    title({'Feed forward neural network',['Nodes in the hidden layer: ' num2str(net.layers{1}.size)]})
    hold on;
    plot3(sonarInput(:,1),sonarInput(:,2),sonarTarget,'bo')
    
    
end