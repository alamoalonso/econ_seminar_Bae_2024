%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code performs the following:
%
% (1) Loads data for the Euro Area or a country of interest. The input data are already de-seasonalized.
% (2) Transforms the data using the transformation codes included in the imput xlsx file
% (3) Imputes missing values and outliers with the EM algorithm (Stock and Watson, 2002; Bańbura and 
%     Modugno, 2014; McCracken and Ng, 2016).
% (4) Imputes data during the Covid perod (2020Q1-2021Q4) for real variables with the EM algorithm,
%     exploting the information included in nominal and financial series.
%
% For additional information on the operations performed by the code and how to choose the relevant 
% parameters, please refer to the ReadME file in the folder.
% ----------------------------------------------------------------------------------------------------------
% OUTPUT:
%        Excel file (.xlsx) with transformed data, outliers and/or Covid period imputed.
% ----------------------------------------------------------------------------------------------------------
% SEE ALSO: EAtransform, remove_outliers, EMimputation
% ----------------------------------------------------------------------------------------------------------
% REFERENCES:
% [1] Bai, J. and Ng, S. "Determining the Number of Factors in Approximate Factor Models." 
%     Econometrica, 2002, 70, 191-221                                                  [BN02]
% [2] McCracken, M.W., and Ng, S. "FRED-MD: A monthly database for macroeconomic research."
%     Journal of Business & Economic Statistics, 34(4), 574-589. (2016)                [MN16]
% [3] Stock, J.H., and Watson, M.W. "Macroeconomic forecasting using diffusion indexes."
%     Journal of Business & Economic Statistics, 20(2), 147-162. (2002)                [SW02]
% [4] Bańbura, M., and Modugno, M. (2014). Maximum likelihood estimation of factor 
%     models on datasets with arbitrary pattern of missing data. 
%     Journal of applied econometrics, 29(1), 133-160.                                 [BM14]
% -----------------------------------------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona@gmail.com)
% Last update: 27/Nov/2023
% Version: MATLAB 2023b
% Required Toolboxes: /
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
parts = strsplit(pwd, filesep); parent_path = strjoin(parts(1:end-1), filesep);
addpath(genpath(parent_path));

if ispc; sep = '\'; else; sep = '/'; end                                                                                 % define separator (for paths)

[sheet,IDf,trans,idt,method,q,opts.algo] = settings(); 

if strcmp(sheet,'all'); country = {'AT','BE','DE','EA','EL','ES','FR','IE','IT','NL','PT'}; 
else; country = {sheet}; 
end

saveto = join([pwd,sep,"data_TR",string(idt)],''); [~,~,~] = mkdir(saveto);
for cc=1:numel(country)
    %%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 1: Load data %%
    %%%%%%%%%%%%%%%%%%%%%%%
    sheet = country{cc};    % select country sheet

    dataset = readtable(join([sheet,'data','.xlsx'],''),'Sheet','data','VariableNamingRule','preserve');                 % upload data
    info = readtable(join([sheet,'data','.xlsx'],''),'Sheet','info');                                                    % upload info
    spec = table2struct(info, "ToScalar",true);                                                                          % convert to structure
    
    spec.TR = spec.(join(['TR',string(idt)],''));                                                                        % select transformations

    data = table2array(dataset(:,2:end));                                                                                % array of data
    titles = string(dataset(:,2:end).Properties.VariableNames);                                                          % titles
    dates = table2array(dataset(:,1));                                                                                   % dates

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 2: Set frequency %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch IDf
     %--------                                                                                                           %=================================%
     case 'QM'                                                                                                           %  CASE 1: quarterly aggregation
     %--------                                                                                                           %=================================%
        [data,dates] = aggregate(data,dates,spec);                                                                       % aggregate 'M' to 'Q'                                                                      
     %--------                                                                                                           %=================================%
     case 'Q'                                                                                                            %   CASE 2: only quarterly data
     %--------                                                                                                           %=================================%
         data = data(:,strcmp(spec.Frequency,'Q'));                                                                      % subset data, titles, period, info      
         titles = titles(:,strcmp(spec.Frequency,'Q'));                                                                  % ---------------------------------
         spec.TR = spec.TR(strcmp(spec.Frequency,'Q'),:);                                                                % ---------------------------------
         spec.Class = spec.Class(strcmp(spec.Frequency,'Q'),:);                                                          % ---------------------------------
         spec.Name = spec.Name(strcmp(spec.Frequency,'Q'),:);                                                            % ---------------------------------
         spec.Aggregation = spec.Aggregation(strcmp(spec.Frequency,'Q'),:);                                              % ---------------------------------
         spec.Frequency = spec.Frequency(strcmp(spec.Frequency,'Q'),:);                                                  % ---------------------------------
         [data,dates] = aggregate(data,dates,spec);                                                                      % aggregate 'M' to 'Q'                                                                      
     %--------                                                                                                           %=================================%
     case 'M'                                                                                                            %    CASE 3: only monthly data
     %--------                                                                                                           %=================================%
         data = data(:,strcmp(spec.Frequency,'M'));                                                                      % subset data, titles, period, info      
         titles = titles(:,strcmp(spec.Frequency,'M'));                                                                  % ---------------------------------
         spec.TR = spec.TR(strcmp(spec.Frequency,'M'),:);                                                                % ---------------------------------
         spec.Class = spec.Class(strcmp(spec.Frequency,'M'),:);                                                          % ---------------------------------
         spec.Name = spec.Name(strcmp(spec.Frequency,'M'),:);                                                            % ---------------------------------
         spec.Aggregation = spec.Aggregation(strcmp(spec.Frequency,'M'),:);                                              % ---------------------------------
         spec.Frequency = spec.Frequency(strcmp(spec.Frequency,'M'),:); 
     %--------                                                                                                           %=================================%
     otherwise; spec.agg = 0;                                                                                            %    CASE 4: MIXED FREQUENCIES    %
     %--------                                                                                                           %=================================%              
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 3: Transform data %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spec.trans = trans;                                                                                                  % select transformations
    xt = EAtransform(data,spec);                                                                                         % transform data

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 4: Impute outliers/missing values/ragged edges %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    opts.maxiter = 1000; opts.thresh = 10^-5;                                                                            % select EM algorithm tuning parameters
    T19 = find(dates=='01-Oct-2019'); T21 = find(dates=='01-Oct-2021');                                                  % locate 2019Q4 and 2021Q4
    opts.T19 = T19;                                                                                                      %=================================%      
    switch method
        %=====                                                                                                           %===============================================%    
        case -1                                                                                                          %      CASE -1: no imputation whatsoever
        %-----                                                                                                           %===============================================%    
        Xc = xt;                                                                                                         % set transformed data
        IDi = 'NA';                                                                                                      % no label
        %=====                                                                                                           %===============================================%    
        case 0                                                                                                           %       CASE 0: impute only ragged edges
        %-----                                                                                                           %===============================================%    
          opts.out = 0;                                                                                                  % no outlier imputation
          Xc = EMimputation(xt,q,opts);                                                                                  % EM algorithm w/out outlier imputation
          IDi = '';                                                                                                      % no label                                         
        %-----                                                                                                           %===============================================%    
        case 1                                                                                                           %      CASE 1: standard outlier imputation
        %-----                                                                                                           %===============================================%    
            opts.out = 1;                                                                                                % set outlier imputation
            [Xc,pc] = EMimputation(xt,q,opts);                                                                           % EM algorithm w/ outlier imputation
            IDi = 'OA';                                                                                                  % id imputation (for saving)
        %-----                                                                                                           %===============================================%    
        case 2                                                                                                           %   CASE 2: impute real variables during Covid
        %-----                                                                                                           %===============================================%    
            loc = find(strcmp(spec.Class,'R'));                                                                          % locate real variables
            Xnan = xt; Xnan(T19+1:T21,loc) = nan;                                                                        % NaN-out Covid period
            [Xc,pc] = EMimputation(Xnan,q,opts);                                                                         % impute outliers w/ [MN16]
            IDi = 'COV';                                                                                                 % id imputation (for saving)
    end

    TT= array2timetable(Xc,'RowTimes',dates,'VariableNames',spec.Name);                                                  % timetable of data
    if ~isempty(IDi); sIDi = join(['_',IDi],''); else; sIDi = IDi; end                                                   % define label: imputation
    nsave = join([saveto,sep,sheet,'data',IDf,sIDi,'_TR',string(idt),'.xlsx'],'');                                       % file name
    writetimetable(TT,nsave,'Sheet',sheet);                                                                              % ---------------------

end

%%%%%%%%%%%%%%%%%
%%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------------------------------------------
function [sheet,IDf,trans,idt,method,q,algo] = settings()
% =========================================================================
% Pop-up box to choose model settings.
% ---------------------------------------------------------------------------------
% SYNTAX: settings()
% ----------------------------------------------------------------------------------
% OUTPUT:
%        - sheet: identifier of the country
%        - IDf: frequency
%        - trans: transformation method
%        - method: imputation method
%        - q = number of factors used for imputation
% ----------------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona2@unibo.it)
% Last update: 01/Dec/2023
% Version: MATLAB 2023b
% Required Toolboxes: /
% ----------------------------------------------------------------------------------
% SEE ALSO:
% -----------------------------------------------------------------------------------

answer = questdlg('Use default options?','Setup','Yes','No','Yes'); % dialogue box

switch answer
  %---------                                                                                                             %================================%    
  case 'Yes'                                                                                                             % 1.       DEFAULT OPTIONS       %    
  %---------                                                                                                             %================================%    
        sheet = 'all';                                                                                                   % COUNTRIES: all
        IDf = 'QM';                                                                                                      % FREQUENCY: Q, aggregate M
        trans = 'light';                                                                                                 % TRANSFORMATION: light
        method = 0;                                                                                                      % IMPUTATION: only ragged edges
        q = 99;                                                                                                          % N. FACTORS: select
        algo = 'SW';                                                                                                     % ALGORITHM: Stock and Watson 
        %-----------
        disp(join(['Country: ' sheet '; ' 'Frequency: ' IDf '; ' 'trans: ' trans '; '...                                 % display choices
             'method: ' string(method) '; ' 'q: ' string(q) '; '...                                                      % ---------------
             'algorithm: ' string(algo) '.'],''))                                                                        % ---------------
        %-----------
  case 'No'
        sheet = inputdlg('Select country','Country',[1 20],{'all'});                                                     % COUNTRIES: select
        sheet = char(sheet);                                                                                             % get selected country
        if isempty(sheet) 
           disp('No country specified, downloading data for all countries');                                             % if empty, set to default
           sheet = 'all';                                                                                                % default: all countries
        end

        IDf = listdlg('PromptString',{'Select frequency'},'SelectionMode','single',...                                   % FREQUENCY: select
              'ListString',{'QM', 'Q', 'M','MF'},'InitialValue',1);                                                      % -----------------
        IDf = char(IDf);                                                                                                 % get selected frequency
        
        trans = questdlg('Select transformation','Transformation','light','heavy','BLT','light');                        % TRANSFORMATION: select

        imp = questdlg('Impute missing values/outliers/Covid period?','Imputation','Yes','No','Yes');                    % IMPUTATION: select 
        imp = char(imp);                                                                                                 % get selected choice

        if strcmp(imp,'Yes')                                                                                             % if imputation selected:
           method = inputdlg('Select imputation method','Imputation',[1 30],{'0'});                                      % IMPUTATION METHOD: select
           method = str2double(method);                                                                                  % get selected method
           
           q = inputdlg('Select number of factors','Factors',[1 25],{'99'});                                             % N. OF FACTORS: select
           q = str2double(q);                                                                                            % get selected number of factors
           
           if isempty(q)                                                                                                 % if no number specified: 
              disp('Number of factors unspecified, setting to default (99)'); q = 99;                                    % select factors (Bai&Ng 2002)
           end

           algo = questdlg('Select imputation algorithm','Imputation Algorithm','SW','BM','SW');                         % ALGORITHM: select

        else                                                                                                             % if no imputation selected...
           method = -1; q = 99; algo = '';                                                                               % ...return data with NaNs 
        end
        %-----------
        disp(join(['Country: ' sheet '; ' 'Frequency: ' IDf '; ' 'trans: ' trans '; '...                                 % display choices
             'method: ' string(method) '; ' 'q: ' string(q) '; '...                                                      % ---------------
             'algorithm: ' string(algo) '.'],''))                                                                        % ---------------
        %-----------
 end

 if strcmp(trans,'light'); idt = 2;                                                                                      % all rates in levels
 elseif strcmp(trans,'heavy'); idt = 1;                                                                                  % statistical transformations
 elseif strcmp(trans,'BLT'); idt = 3;                                                                                    % transformations in [BLT25]
 end

 if strcmp(IDf,'MF')&method>-1                                                                                           % if MF + imputation...
    warning('Imputation of missing values not available for mixed-frequency data, setting <method> to -1')               % ...throw warning
    method = -1; 
 end

end


%-----------------------------------------------------------------------------------------------------------------------

function [dataQ,datesQ,datasetQ] = aggregate(data,dates,spec)
% ----------------------------------------------------------------------------------
% Aggregate a mixed-frequency monthly-quarterly dataset to quarterly levels. 
% ----------------------------------------------------------------------------------
% SYNTAX: aggregate(data,dates,spec)
% where:  
%      - data = matrix of data with NaNs (TxN)  
%      - dates = vector of dates (Tx1 datetime)
%      - spec = structure containing:
%              Frequency: vector of strings identifying quarterly ('Q') or
%                         monthly ('M') variables (Nx1)
%              Aggregation: array of identifiers for aggregation method (Nx1)
%                          = 1 if aggregation by mean (stock variables)
%                          = 2 if aggregation by sum  (flow variables) 
%              Name: array of variable names (Nx1)
%                    (optional)
% ----------------------------------------------------------------------------------
% OUTPUT: 
%      - dataQ = matrix of quarterly data (TqxN)
%      - datesQ =  vector of quarterly dates (Tqx1)
%      - datasetQ = table with headers and dates as row index (TqxN)
% ----------------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona2@unibo.it)
% Date: 24/Aug/2023
% Version: MATLAB 2021b
% Required Toolboxes: /
% -----------------------------------------------------------------------------------

temp = data;                                                                                                             % store de-aggregated data
tm = median(dates);                                                                                                      % meadian date 

if numel(spec.Frequency)~=size(data,2)||numel(spec.Aggregation)~=size(data,2) % check size mismatch...
    error('Identifiers must have dimension N'); % ... and throw error
end

for i = find(strcmp(spec.Frequency,'M'))'                                                                                % iterate over monthly variables
  
    idxNaN = isnan(data(:,i));                                                                                           % locate NaNs
    t0 = min(dates(~idxNaN)); t1 = max(dates(~idxNaN));                                                                  % min/max non-NaN date
    if t0==t1; if t0<tm; t1=[]; else; t0=[]; end; end                                                                    % if min=max, chack whether is > or < of median date
    
    if ~isempty(t0)                                                                                                      % min non-NaN date
        tp = dateshift(t0,'start','month','next');                                                                       % month + 1
        tp1 = dateshift(tp,'start','month','next');                                                                      % month + 2
        if (quarter(tp)~=quarter(t0))||(quarter(tp1)~=quarter(t0))                                                       % if one value is missing, set as NaN 
           temp(dates==t0,i) = nan;                                                                                      % ---> need three full months to build quarterly value
        end                                                                                                              %-----------------------------------------------------
    end

    if ~isempty(t1)                                                                                                      % max non-NaN date
        tp = dateshift(t1,'start','month','previous');                                                                   % month - 1
        tp1 = dateshift(tp,'start','month','previous');                                                                  % month - 2
        if (quarter(tp)~=quarter(t1))||(quarter(tp1)~=quarter(t1))                                                       % if one value is missing, set as NaN 
           temp(dates==t1,i) = nan;                                                                                      % ---> need three full months to build quarterly value
        end                                                                                                              %-----------------------------------------------------
    end
end

data = temp;                                                                                                             % update date               

%---------------------%
% AGGREGATE QUARTERLY %
%---------------------%
if ~isfield(spec,'Name'); spec.Name = cellstr(string((1:size(data,2))')); end                                            % if no names provided, use numbers as column names
TT = array2timetable(data,'RowTimes',dates,'VariableNames',spec.Name);                                                   % build table
if any(strcmp(spec.Frequency,'M'))
   id_flow = find(strcmp(spec.Frequency,'M')&spec.Aggregation==2);                                                       % identify flow variables
else; id_flow = [];
end
TT2 = retime(TT,'quarterly',@(x) mean(x,'omitnan'));                                                                     % stock variables: aggregate by taking means
TT2(:,id_flow) = retime(TT(:,id_flow),'quarterly',@(x) sum(x));                                                          % flow variables: aggregate by taking sum

datasetQ = timetable2table(TT2); dataQ = table2array(datasetQ(:,2:end));                                                 % dataset (table) and data (array)
datesQ = table2array(datasetQ(:,1));                                                                                     % dates
end

%-----------------------------------------------------------------------------------------------------------------------
function Xt = EAtransform(X,spec,c)
% =========================================================================
% Transform dataset of variables given an indicator of transformations (TR). 
% For each variable x:
%
%      - TR = 1 ---> c*log(x)
%      - TR = 2 ---> \Delta log(c*x)
%      - TR = 3 ---> \Delta(\Delta log(c*x))
%      - TR = 4 ---> x (no transformation)
%      - TR = 5 ---> \Delta x
%      - TR = 6 ---> \Delta(\Delta x)
%
% whenever TR = 1.5, 2.5, 4.5, the transformation depends on the
% methodology: if light transformations, take floor(TR), if heavy
% transformations, take ceil(TR). The function allows for the presence
% of missing values due to, e.g., mixed frequencies.
% ---------------------------------------------------------------------------------
% SYNTAX: EAtransform(X,spec,c)
% where:
%      - X = matrix of data (TxN)
%      - spec = structure containing
%               TR: vector of identifiers for transformations (Nx1)
%                   (admissible values: 0 to 5)
%               trans: 1 for light transformations, 2 for heavy transformations
%               Frequency: vector of strings identifying quarterly ('Q') or
%                          monthly ('M') variables (Nx1)
%               agg: 1 if data are aggregated, 0 if mixed frequencies to account
%                    (optional, default is 1)
%    - c = multiplicative factor for differences: if 100, takes %growth rates.
%          (optional, default is 1)
% ----------------------------------------------------------------------------------
% OUTPUT:
%        - Xt =  matrix of transformed data (TxN)
% ----------------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona2@unibo.it)
% Last update: 28/Aug/2023
% Version: MATLAB 2023b
% Required Toolboxes: /
% ----------------------------------------------------------------------------------
% SEE ALSO: mdiff
% -----------------------------------------------------------------------------------
if nargin < 3; c = 1; end % default: do not scale differences
if ~isfield(spec,'agg'); spec.agg = 1; end                                                                               % default: aggregate

if size(spec.TR,1)~=size(X,2)                                                                                            % check size
    if size(spec.TR,1)==size(X,1); X = X';                                                                               % ----------
    else; error('The number of series in x must coincide with the size of TR'); end                                      % ----------
end                                                                                                                      

if sum(any((X<0)&(spec.TR'==1|spec.TR'==2|spec.TR'==3)))>0                                                               % check whether logs are admissible
    loc = find(any((X<0)&(spec.TR'==1|spec.TR'==2|spec.TR'==3)));                           % identify unadmissible tranformations
    warning(join(['Variables in positions ' + join(string(loc),', ') + ' contain negative values \n'...
                  'logs are not admitted.','Transforming without logs'],''))                           % display error
    if spec.TR(loc)==1; spec.TR(loc) = 4;
    elseif spec.TR(loc)==2; spec.TR(loc) = 5;
    elseif spec.TR(loc)==3; spec.TR(loc) = 6;
    end
end                                                                                                                      %

Xt = NaN(size(X));                                                                                                       % initialize output matrix

%-----------------%
% TRANSFORMATIONS %
%-----------------%
if ~any(spec.agg==[0,1]); spec.agg = 1; disp('Invalid value for <agg>, setting to 1'); end
%-----------------------------------------------------------------------------------------------------------------       %===============================
if spec.agg==1                                                                                                           % BALANCED PANEL
    %-----------------------------------------------------                                                               %-------------------------------
    Xt(:,spec.TR==1) = c*log(X(:,spec.TR==1));                                                                           % logs
    Xt(2:end,spec.TR==2) = diff(c*log(X(:,spec.TR==2)));                                                                 % log-differences
    Xt(3:end,spec.TR==3) = diff(c*log(X(:,spec.TR==3)),2);                                                               % second log-differences
    Xt(:,spec.TR==4) = X(:,spec.TR==4);                                                                                  % no transformation
    Xt(2:end,spec.TR==5) = diff(X(:,spec.TR==5));                                                                        % first-differences
    Xt(3:end,spec.TR==6) = diff(X(:,spec.TR==6),2);                                                                      % second-differences
    %-----------------------------------------------------                                                               %--------------------------------
elseif spec.agg==0                                                                                                       % MIXED-FREQUENCIES PANEL
    %-------------------------------------------                                                                         %--------------------------------
    locQ = strcmp(spec.Frequency,'Q'); locM = strcmp(spec.Frequency,'M');   % locate Q and M variables

    Xt(:,spec.TR==1) = c*log(X(:,spec.TR==1));                                                                           % logs
    Xt(2:end,locM&spec.TR==2) = diff(c*log(X(:,locM&spec.TR==2)));           % log-differences: M vars
    Xt(:,locQ&spec.TR==2) = mdiff(c*log(X(:,locQ&spec.TR==2)),3,1);          % log-differences: Q vars
    Xt(3:end,locM&spec.TR==3) = diff(c*log(X(:,locM&spec.TR==3)),2);         % second log-differences: M vars
    Xt(:,locQ&spec.TR==3) = mdiff(c*log(X(:,locQ&spec.TR==3)),3,2);          % second log-differences: Q vars
    Xt(:,spec.TR==4) = X(:,spec.TR==4);                                                                                  % no transformation
    Xt(2:end,locM&spec.TR==5) = diff(X(:,locM&spec.TR==5));                  % first-differences: M vars
    Xt(:,locQ&spec.TR==5) = mdiff(X(:,locQ&spec.TR==5),3,1);                 % first-differences: Q vars
    Xt(3:end,locM&spec.TR==6) = diff(X(:,locM&spec.TR==6),2);                % second-differences: M vars
    Xt(:,locQ&spec.TR==6) = mdiff(X(:,locQ&spec.TR==6),3,2);                 % second-differences: Q vars
    %-------------------------------------------                                                                         %--------------------------------
end                                                                                                                      %
%-----------------------------------------------                                                                         %================================
end

%-----------------------------------------------------------------------------------------------------------------------
function pc = princfact(X,q,method,st)
% ==========================================================================
% Function to estimate factors via principal components. The functions
% allows for the presence of NaNs, which are replaced by the sample mean of
% non-missing observations
% --------------------------------------------------------------------------
% SYNTAX: princfact(X,q,method,st)
% where:
%      - X = matrix of data (TxN)
%      - q = number of factors to extract
%      - method = normalization of the factors
%               0: no normalization
%               1: scale by eigvals
%               2: scale by N (default)
%               3: estimate using TxT var-cov matrix
%      - st = 1 to standardize the panel, 0 otherwise (default: st=0)
% ---------------------------------------------------------------------------
% OUTPUT: pc (structure)
% where:
%      - F = matrix of estimated factors (Txq)
%      - C = matrix of estimated loadings (Nxq)
%      - chi = estimated common component (TxN)
%      - R = estimated idosyncratic variances (TxN,diagonal)
%      - d = q-largest eigenvalues (qx1)
% ---------------------------------------------------------------------------
% REFERENCES:
%
% Barigozzi, M. "On Estimation and Inference of Large Approximate Dynamic
%                Factor Models via the Principal Component Analysis."
%                arXiv preprint arXiv:2211.01921 (2022).                 [B23]
% ----------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona2@unibo.it)
% Last update: 24/Aug/2023
% Version: MATLAB 2021b
% Required Toolboxes: /
% ----------------------------------------------------------------------------
if nargin<3; method=2; end                                                                                               % default method: 2
if nargin<4; st=0; end                                                                                                   % default std: none

[T,N] = size(X);                                                                                                         % size of the panel

if any(isnan(X))                                                                                                         %------------
    X(isnan(X)) = mean(X(~isnan(X)));                                                                                    % replace NaNs with the unconditional mean of
    disp('Detected NaNs in original data matrix: replaced with unconditional mean of non-NaNs values')                   % each series and display warning
end                                                                                                                      %------------

if st==1                                                                                                                 %------------------
    X = (X - repmat(mean(X),1,size(X,1)))./repmat(std(X),1,size(X,1));                                                   % standardize panel
end                                                                                                                      %------------------

if ~any(method==[0,1,2,3]); disp('The value for <method> is outside the bounds, setting <method> to 2'); end             % if method uncorrectly specified, set to 2

%-----------------%
% Factors via PCA %
%-----------------%
[v,d] = eigs(cov(X),q,'LM');                                                                                             % extract q-largest eigvalues and associated eigvectors
if(sum(v) < 0); v = -v; end                                                                                              % flip sign (just for better interpretability)
%----------------------------------                                                                                      %--------------------------------------------
if method==0                                                                                                             % METHOD 1: C'C = I_q; cov(F) diagonal, with eigenvalues as variances --->
    f = X*v;                                                                                                             % factors
    C = v;                                                                                                               % loadings
    %------------------------------                                                                                      %----------
elseif method==1                                                                                                         % METHOD 2: C'C diagonal, with eigenvalues as variances; cov(F) = I_q ---> see [B23] Section 2.1
    f = X*v*diag(1./diag(sqrt(d)));                                                                                      % factors
    C = v*sqrt(d);                                                                                                       % loadings
    %------------------------------                                                                                      %----------
elseif method==2                                                                                                         % METHOD 3: C'C diagonal, with variances equal to N; cov(F) = eigenvalues scaled by N ---> see [B23] Section 2.3
    f = X*v/sqrt(N);                                                                                                     % factors
    C = v*sqrt(N);                                                                                                       % loadings
    %------------------------------                                                                                      %-----------
elseif method==3                                                                                                         % METHOD 4: cov(C) diagonal with eigenvalues as variances; cov(F) =  I_q ---> see [B23] Section 2.2
    [v,d] = eigs(cov(X'),q,'LM');                                                                                        % use TxT var-cov matrix
    if(sum(v) < 0); v = -v; end                                                                                          % flip sign (just for better interpretability)
    f = v*sqrt(T);                                                                                                       % factors
    C = (f'*X)'/T;                                                                                                       % loadings
end                                                                                                                      %------------
%----------------------------------                                                                                      %--------------------------------------------
chi = f*C';                                                                                                              % common component
R = diag(diag(cov(X-chi)));                                                                                              % variance idiosyncratic component ---> diagonal (approximate factor model)

pc.F = f; pc.C = C; pc.chi = chi; pc.R = R; pc.d = diag(d);  % build results
end

%-----------------------------------------------------------------------------------------------------------------------

function [X,out,n] = remove_outliers(X,c)
% =========================================================================
% Function to locate outliers and replace them with NaNs
%
% An outlier is identified for each series whenever |abs(X-m_X)>c*iqr_X|,
% where:
%      - m_X = median of X
%      - iqr_X = interquartile range of X
%      - c = multiplicative constant
% -------------------------------------------------------------------------
% SYNTAX: remove_outliers(X,c)
% where:
%      - X = matrix of data with NaNs (TxN)
%      - c = scale of identification: the higher is c, the lower the number
%            of outliers (default: c=10)
% -------------------------------------------------------------------------
% OUTPUT:
%      - Y = matrix of data with outliers as NaNs (TxN)
%      - out = position of outlier (TxN, logical)
%      - n = number of outliers
% -------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona2@unibo.it)
% Last update: 29/Aug/2023
% Version: MATLAB 2021b
% Required Toolboxes: /
% -------------------------------------------------------------------------
if nargin < 2; c = 10; end                                                                                               % default threshold: 10

mX = median(X,1,'omitnan');                                                                                              % median of each series
iqrX = iqr(X);                                                                                                           % interquartile range of each series
out = abs(X-repmat(mX,size(X,1),1)) > c*repmat(iqrX,size(X,1),1);                                                        % identify outliers  % locate col/row of outliers
n = sum(out,1);                                                                                                          % number of outliers

loc = sum(out,1)/size(X,1)>0.2; out(:,loc) = false;    % control for series with many zeroes ---> they introduce "artificial" outliers

X(out) = NaN;                                                                                                            % replace outliers with NaN

end

%-----------------------------------------------------------------------------------------------------------------------

function Y = mdiff(X,s,k)
% =========================================================================
% Function to perform k-th differences with a mixed-frequency panel of data
% (along axis 1)
% -------------------------------------------------------------------------
% SYNTAX: mdiff(X,s,k)
% where:
%      - X = matrix of data with NaNs (TxN)
%      - s = period length (s=3 for quarters, s=12 for years)
%      - k = order of differencing (default: k=1)
% -------------------------------------------------------------------------
% OUTPUT:
%      - Y = matrix of differenced data (TxN)
% -------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona2@unibo.it)
% Last update: 28/Aug/2023
% Version: MATLAB 2023b
% Required Toolboxes: /
% -------------------------------------------------------------------------

if nargin < 2; s = 1; end                                                                                                % if no s and k, boils down to diff(X) with initial NaNs
if nargin < 3; k = 1; end                                                                                                % if no k, take 1-st differences

Y = nan(size(X));                                                                                                        % initialize output matrix
dX = X(s+1:end,:) - X(1:end-s,:);                                                                                        % 1-st differences between periods s

%--------------------------------------                                                                                  %=========================
if k>1                                                                                                                   % Higher-order differences
    for i=1:k-1                                                                                                          %-----------------
        dX = dX(s+1:end,:) - dX(1:end-s,:);                                                                              % k-th differences
    end                                                                                                                  %-----------------
end                                                                                                                      %
%--------------------------------------                                                                                  %==========================
Y(s*k+1:end,:) = dX;                                                                                                     % build output matrix

end

%-----------------------------------------------------------------------------------------------------------------------

function q = BaiNg(X,qmax,method)
% =========================================================================
% Function to select the number of factors in a T x N panel of series X 
% according to one of the information criteria by Bai and Ng (2002), ECMA.
% 
% The model for X_{t} N-dimensional is given by:
%   
%       X_{t} = L F_{t} + e_{t}
%
% where F_{t} is the q-dimensional vector of factors. 
%
% The function allows for the presence of missing values, which will be 
% replaced by the unconditional mean of the corresponding series.
% -------------------------------------------------------------------------
% SYNTAX: BaiNg(X,qmax,method)
% where:
%      - X = matrix of data (T x N)
%      - qmax = maximum number of factors to test 
%               (optional, default is sqrt(N))
%      - method = log-criterion 1, 2 or 3 (optional, default is 2)
% -------------------------------------------------------------------------
% OUTPUT:
%      - q = number of factor selected using <method>
% -------------------------------------------------------------------------
% REFERENCES:
%      1. Bai and Ng (2002). Determining the Number of Factors in 
%         Approximate Factor Models. Econometrica.                   [BN02]
% -------------------------------------------------------------------------
% COMMENTS: Adapted from the original function of the authors. All mistakes
%           are mine.
% -------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona@gmail.com)
% Last update: 19/06/2024
% Version: MATLAB 2024b
% Required Toolboxes: /
% -------------------------------------------------------------------------
[T,N] = size(X);                                                                                                         % size of the panel

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (0) Options and data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2; qmax = sqrt(N); end                                                                                       % default max. q: sqrt(N)
if nargin < 3; method = 2; end                                                                                           % default criterion: 2

if any(isnan(X))                                                                                                         
    X(isnan(X)) = mean(X(~isnan(X)));                                                                                    % replace NaNs with the unconditional mean of
    disp('Detected NaNs in original data matrix: replaced with unconditional mean of non-NaNs values')                   % each series and display warning
end  

X = (X - ones(T,1)*(sum(X,1)/T))./std(X);                                                                                % standardize data

CT = zeros(1,qmax);                                                                                                      % initialize criterion

%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) Select penalty %%
%%%%%%%%%%%%%%%%%%%%%%%%
switch method                                                                                                           
    %------                                                                                                              %========================   
    case 1                                                                                                               %         CASE 1
    %------                                                                                                              %========================   
        CT(1,:) = log(N*T/(N+T))*(1:1:qmax)*(N+T)/(N*T);                                                                 % penalty PC_p1 in [BN02]
    %------                                                                                                              %========================   
    case 2                                                                                                               %         CASE 2
    %------                                                                                                              %========================   
        CT(1,:) = ((N+T)/(N*T))*log(min(N,T))*(1:1:qmax);                                                                % penalty PC_p2 in [BN02]
    %------                                                                                                              %========================   
    case 3                                                                                                               %         CASE 3
    %------                                                                                                              %========================   
        CT(1,:) = (1:1:qmax)*log(min(N,T))/min(N,T);                                                                     % penalty PC_p3 in [BN02]
end                                                                                                                      %========================   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) Principal components %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if T < N
    [egvc,~] = eigs(cov(X',1),qmax,'LM');                                                                                % if T < N, use T x T var-cov matrix

    Fhat = sqrt(T)*egvc;                                                                                                 % factors
    Lhat = X'*Fhat/T;                                                                                                    % loadings
else
    [egvc,~] = eigs(cov(X,1),qmax,'LM');                                                                                 % if T > N, use N x N var-cov matrix

    Lhat = sqrt(N)*egvc;                                                                                                 % factors
    Fhat = X*Lhat/N;                                                                                                     % loadings
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) Information Criterion %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma = zeros(1,qmax+1);                                                                                                 % initialize
IC1 = zeros(1,qmax+1);                                                                                                   % ----------

for qq=qmax:-1:1                                                                                                         % loop over n. of factors until qmax

    chat = Fhat(:,1:qq)*Lhat(:,1:qq)';                                                                                   % common component using qq factors 
    ehat = X - chat;                                                                                                     % idiosyncratic component using qq factors

    Sigma(qq) = mean(sum(ehat.*ehat/T));                                                                                 % sum of squared residuals squared by NT

    IC1(qq) = log(Sigma(qq)) + CT(qq);                                                                                   % information criterion for qq factors ([BN02] eq. 9)

end

Sigma(qmax+1) = mean(sum(X.*X/T));                                                                                       % sum of squared residuals with no factors
IC1(qmax+1) = log(Sigma(qmax+1));                                                                                        % information criterion with no factors

%%%%%%%%%%%%%%%%%%%%%%%
%% (4) N. of factors %%
%%%%%%%%%%%%%%%%%%%%%%%
[~,q] = min(IC1);                                                                                                        % n. of factors minimizes IC
q = q*(q<=qmax);                                                                                                         % if q > qmax then q = 0


end

%-----------------------------------------------------------------------------------------------------------------------

function KF = kalman(Y,pars)
% ----------------------------------------------------------------
% Apply Kalman Filter and Smoother to a panel of time series with
% potential missing data. The model takes the form:
%  
% Y_{t} = beta_1 + beta_2*t + C Z_t + e_{t},     e_{t} ~ N(0,R),
% Z_{t} = mu + A Z_{t-1} + u_{t},                u_{t} ~ N(0,Q)
% 
% with: 
%      Z_{t}|y_{t-1} ~ N(Z_{t|t-1},P_{t|t-1})
%      Y_{t}|Y^{t-1} ~ N(Y_{t|t-1},C*P_{t|t-1}*C'+R)
% ----------------------------------------------------------------
% SYNTAX: KF = KF(Y,pars)
% where:  
%       Y = matrix of data with NaNs (TxN)  
%       pars = structure containing:
%              C: matrix of loadings (Nxns)
%              R: var-cov matrix of idio (NxN)   
%              A: transition matrix (nsxns)
%              Q: var-cov matrix of latent errors (nsxns)
%              Z00: initial vector of states (Tx1)
%              P00: initial var-cov states (nsxnsxT)
%------------------------------------------------------------------
% OUTPUT: KF (structure)
% with:
%     KF:
%       Zttm: predicted states Z_{t|t-1} (ns x T)
%       Pttm: predicted var-cov states P_{t|t-1} (ns x ns x T)
%       Ztt: filtered states Z_{t|t} (ns x (T+1))
%       Ptt: filtered var-cov states P_{t|t} (ns x ns x (T+1))
%       Kt: Kalman gain (with 0s as missings) (ns x N x T)
%       vt: prediction error (with 0s as missings) (N x T)
%       loglik: log-likelihood from Kalman filter
%       ZtT: smoothed states Z_{t|T} (ns x T)
%       PtT: smoothed var-cov stated P_{t|T} (ns x ns x T)   
%       PtTm: smoother autocovariance (ns x ns)
% ------------------------------------------------------------------
% REFERENCES:
%      ---------------------------------------------------------------------------
%  [1] Harvey, Andrew C. "Forecasting, Structural Time Series Models and
%      the Kalman Filter". Cambridge Books, Cambridge University Press.      [H90]
% ---------------------------------------------------------------------------------

A = pars.A; C = pars.C; R = pars.R; Q = pars.Q; 

if isfield(pars,'Z00'); Z0 = pars.Z00; else; Z0 = zeros(size(C,2),1); end                                                % if no initial state, initialize with 0
if isfield(pars,'P00'); P0 = pars.P00; else; P0 = 10*eye(size(C,2)); end                                                 % if no initial state var-cov, initialize with large variance

if isfield(pars,'mu'); mu = pars.mu; else; mu=zeros(size(C,2),1); end                                                    % if no constant in state equation, set to 0
if isfield(pars,'beta'); beta = pars.beta; else; beta=zeros(1,2); end                                                    % if no deterministic term in observation equation, set to 0

[~,ns] = size(C);                                                                                                        % number of states

Y = Y'; [N,T] = size(Y);                                                                                                 % size of the panel

%%%%%%%%%%%%%%%%%%%%%%%%
%% (0) Initialization %%
%%%%%%%%%%%%%%%%%%%%%%%%
KF.Zttm = nan(ns,T); KF.Pttm = nan(ns,ns,T);                                                                             % Z_{t|t-1} and P_{t|t-1}
KF.Ztt = nan(ns,T+1); KF.Ptt = nan(ns,ns,T+1);                                                                           % Z_{t|t} and % P_{t|t} 
KF.Kt = zeros(ns,N,T);  KF.vt = zeros(N,T);                                                                              % Kalman gain and prediction error
KF.loglik = 0;                                                                                                           % log-likelihood
KF.Ztt(:,1) = Z0; KF.Ptt(:,:,1) = P0;                                                                                    % Z_{0} and V_{0} 
mu = [zeros(ns,1), repmat(mu,1,T-1)];                                                                                    % set initial mu to 0 

%%%%%%%%%%%%%%%%%%%%%%%
%% (1) Kalman Filter %%
%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:T
     %-----------------%
     % Prediction Step %  
     %-----------------%                                                                                                 % ---------------                                                              
     KF.Zttm(:,t) = mu(:,t) + A*KF.Ztt(:,t);                                                                             % Z_{t|t-1}
     KF.Pttm(:,:,t) = A*KF.Ptt(:,:,t)*A' + Q;                                                                            % P_{t|t-1}
     KF.Pttm(:,:,t) = 0.5*(KF.Pttm(:,:,t)+KF.Pttm(:,:,t)');                                                              % ensure symmetry
     
     %-------------%                                                                                                     
     % Update Step %  
     %-------------%                                                                                                    
     idx = ~isnan(Y(:,t));                                                                                               % identify series without NaN
     yt = Y(idx,t); Ct = C(idx,:); Rt = R(idx,idx); if(sum(beta(:)))~=0; bt = beta(idx,:); else; bt=beta; end            % remove NaNs from parameters in the obs. equation

     %----------------------------------------------------------------------------                                       % ---------------      
     if isempty(yt)                                                                                                      % if no data, no update ---> use prediction values
     %----------------------------------------------------------------------------                                       % --------------- 
          KF.Ztt(:,t+1) = KF.Zttm(:,t);                                                                                  % Z_{t|t}
          KF.Ptt(:,:,t+1) = KF.Pttm(:,:,t);                                                                              % P_{t|t}
     %----------------------------------------------------------------------------                                       % --------------- 
     else                                                                                                                % if data, update
     %----------------------------------------------------------------------------                                       % --------------- 
          KF.Kt(:,idx,t) = KF.Pttm(:,:,t)*Ct'/(Ct*KF.Pttm(:,:,t)*Ct' + Rt);                                              % Kalman gain 
          KF.vt(idx,t) = yt - Ct*KF.Zttm(:,t) - bt*[1;t];                                                                % prediction error v_{t|t-1}
          KF.Ztt(:,t+1) = KF.Zttm(:,t) + KF.Kt(:,:,t)*KF.vt(:,t);                                                        % Z_{t|t}
          KF.Ptt(:,:,t+1) = KF.Pttm(:,:,t) - KF.Kt(:,idx,t)*Ct*KF.Pttm(:,:,t);                                           % P_{t|t}
          KF.Ptt(:,:,t+1) = 0.5*(KF.Ptt(:,:,t+1)+KF.Ptt(:,:,t+1)');                                                      % ensure symmetry
          
          KF.loglik = KF.loglik + 0.5*(log(det(inv(Ct*KF.Pttm(:,:,t)*Ct'+Rt))) -...                                      % compute log-likelihood
                      KF.vt(idx,t)'/(Ct*KF.Pttm(:,:,t)*Ct'+Rt)*KF.vt(idx,t));                                            % ----------------------
     %-----------------------------------------------------------------------------                                      
     end
end
KF.Z00 = Z0; KF.P00=P0;                                                                                                  % initial values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  (2) KALMAN SMOOTHER    %%                                                                                                    
% (FIXED-INTERVAL SMOOTHER) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KF.ZtT = zeros(ns,T); KF.PtT = zeros(ns,ns,T); KF.PtTm = zeros(ns,ns,T);                                                 % Z_{t|T}, P_{t|T}, P_{t,t-1|T}
                                                                                                   
KF.ZtT(:,T) = KF.Ztt(:,T+1); KF.PtT(:,:,T) = KF.Ptt(:,:,T+1);                                                            % initial smoothed states and variances
KF.PtTm(:,:,T) = (eye(ns) - KF.Kt(:,:,T)*C)*A*KF.Ptt(:,:,T);                                                             % initial smoothed autocovariance
Pts = zeros(ns,ns,T);                                                                                                    %  
Pts(:,:,T) = KF.Ptt(:,:,T)*A'*pinv(KF.Pttm(:,:,T));                                                                      % initial update matrix

for t=T:-1:2
    KF.ZtT(:,t-1) = KF.Ztt(:,t) + Pts(:,:,t)*(KF.ZtT(:,t) - KF.Zttm(:,t));                                               % Z_{t|T}  
    KF.PtT(:,:,t-1) = KF.Ptt(:,:,t) + Pts(:,:,t)*(KF.PtT(:,:,t) - KF.Pttm(:,:,t))*Pts(:,:,t)';                           % P_{t|T}
    if t>1
       Pts(:,:,t-1) = KF.Ptt(:,:,t-1)*A'*pinv(KF.Pttm(:,:,t-1));                                                         % update matrix 
       KF.PtTm(:,:,t-1) = KF.Ptt(:,:,t)*Pts(:,:,t-1)'+Pts(:,:,t)*(KF.PtTm(:,:,t)-A*KF.Ptt(:,:,t))*Pts(:,:,t-1)';         % P_{t,t-1|T}
    end
end
 KF.ZtT = KF.ZtT';                                                                                                       % smoothed states
end

%-----------------------------------------------------------------------------------------------------------------------

function res = estVAR(Y,p,opts)
% ==========================================================================
% Estimate VAR model. The model writes as:
%
%  Y_{t} = \gamma_{0} + \gamma_{1}t + \gamma_{2}t^{2} + 
%           + A_{1}Y_{t-1} + ... + A_{p}Y_{t-p} + u_{t}
%  
% In compact form:
%
%  Y_{t} = B[Y_{t-1},...,Y_{t-p}] + u_{t}
%  
% --------------------------------------------------------------------------
% SYNTAX: VAR(Y,p,opts)
% where:  
%      - Y = matrix of data (T x N)  
%      - p = number of lags (optional, estimated if not specified)
%      - opts = structure containing:
%              spec: 0 = none; 1 = constant, 2 = constant + trend, 
%                    3 = constant+ trend + quadratic trend.
%                    (optional, default is 0)
%              crit: if number of lags is estimated, choose criterion
%                    1 = BIC; 2 = AIC; 3 = HQC;
% ---------------------------------------------------------------------------
% OUTPUT: res (structure)
% where:
%               k = number of parameters
%               B: full matrix of coefficients (k x N)
%               Bols: full matrix of coefficients: eq-by-eq (k x N)
%               A: autoregressive coefficients (N*p x N)
%               gamma: coefficients deterministic components (nd x 1)
%               J: companion form autoregressive matrix (N*p x N*p)
%               u: residuals (TxN)
%               Q: variance residuals (N x N)
%               Q_ML: vairance residuals (Maximum-Likelihood) (N x N)
%               loglik: estimated loglikelihood 
% ---------------------------------------------------------------------------
% REFERENCES:
%
% [1] Lütkepohl, H. "New introduction to multiple time series analysis. 
%     Springer Science & Business Media", (2005).                       [L05]
% ---------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona2@unibo.it)
% Last update: 29/Sep/2023
% Version: MATLAB 2023b
% Required Toolboxes: /
% ----------------------------------------------------------------------------

if nargin < 3; opts.spec = 0; end; if ~isfield(opts,'spec'); opts.spec = 0; end                                          % default model: no constant

if ~isfield(opts,'crit'); opts.crit = 1; end                                                                             % default criterion (if any): BIC 
if (opts.crit<1)||(opts.crit>3)
    opts.crit = 1;
    warning('Admissible values for <crit> are 1 to 3, setting to 1 (BIC)')
end

X = lagvar(Y,p,2);                                                                                                       % lagged Y_{t}: all/eq-by-eq  
Y = Y(p+1:end,:);                                                                                                        % Y_{t} 
[T,N] = size(Y);                                                                                                         % size of the panel

%--------------------------------------------------                                                                      % set specification                                           
switch opts.spec                                                                                                         %=====================================
 case 1; X = [ones(T,1),X];                                                                                              % 1. constant
 case 2; X = [ones(T,1),(1:1:T)',X];                                                                                     % 2. constant + trend
 case 3; X = [ones(T,1),(1:1:T)',((1:1:T).^2)',X];                                                                       % 3. constant + linear&quadratic trend
end                                                                                                                      %=====================================
%--------------------------------------------------

res.k = size(X,2); nd = res.k-p*N;                                                                                       % n. of variables/deterministic components

res.B = (X'*X)\(X'*Y);                                                                                                   % coefficients 
res.u = Y - X*res.B; res.Q = (res.u'*res.u)/(T-res.k);                                                                   % residuals and variance  
res.A = res.B(nd+1:end,:); res.gamma = res.B(1:nd,:);                                                                    % autoregressive coefficients/deterministic components 

                                                                                                                         % companion form
if p>1                                                                                                                   %===============
   res.J = zeros(N*p,N*p); res.J(1:N,:) = res.A'; res.J(N+1:end,1:N*(p-1)) = eye(N*(p-1));                               % p>1 
else; res.J = res.A';                                                                                                    % p=1
end                                                                                                                      %===============
 
res.AL = nan(N,N,p); j = 1;          
%---------------------------------------                                                                                 % split autoregressive coefficients                       
for i=0:N:N*(p-1)                                                                                                        %=====================================
    res.AL(:,:,j) = res.A(i+1:i+N,:)';                                                                                   % lag i coefficients
    j = j+1;                                                                                                             % ------------------
end                                                                                                                      %=====================================                      
%---------------------------------------
                                                                                                                         % MAXIMUM-LIKELIHOOD ESTIMATION
                                                                                                                         %===============================
res.Q_ML = res.Q*((T-res.k)/T);                                                                                          % var-cov matrix --> closed-form

res.loglik = -(T*N/2)*log(2*pi) - (T/2)*log(det(res.Q_ML));                                                              % log-likelihood
for t=1:T; res.loglik = res.loglik - (1/2)*(res.u(t,:)/res.Q_ML*res.u(t,:)'); end                                        % --------------

end

%-----------------------------------------------------------------------------------------------------------------------

function X = lagvar(Y,p,method,order)
% ====================================================================
% Obtained lagged variable across the first axis 
% --------------------------------------------------------------------
% SYNTAX: lagvar(Y,p,method)
% where:  
%      - Y = matrix of data (TxN)  
%      - p = number of lags to be considered
%      - method = 1: returns lag 0; 2: returns only positive lags
%      - order = 1: lags the vector; 2: lags each variable sequentially 
%                (optional, default is 1)   
% ---------------------------------------------------------------------
% OUTPUT: 
%      - X = matrix of lagged data
%            > method = 1: (T-p) x N*(p+1), containing also lag 0
%            > method = 2: (T-p) x N*p, without lag 0
% ---------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona2@unibo.it)
% Last update: 27/Sep/2023
% Version: MATLAB 2023b
% Required Toolboxes: /
% ---------------------------------------------------------------------

N = size(Y,2);

if nargin<4, order=1; end
                                                                                                                         % LAG MATRIX      
switch method                                                                                                            %================== 
 case 1                                                                                                                  % Case 1: use lag 0
      %----------------------------------------------------                                                              %------------------ 
      if order==1; X = []; for i=0:p; X = cat(2,X,Y(p-i+1:end-i,:)); end                                                               
      else; X = []; for i=1:N; for j=0:p; X = cat(2,X,Y(p-j+1:end-j,i)); end; end
      end
       %----------------------------------------------------                                                             %------------------
  case 2                                                                                                                 % Case 2: no lag 0
      %----------------------------------------------------                                                              %------------------
      if order==1; X = []; for i=1:p; X = cat(2,X,Y(p-i+1:end-i,:)); end
      else; X = []; for i=1:N; for j=1:p; X = cat(2,X,Y(p-j+1:end-j,i)); end; end
      end  
      %----------------------------------------------------                                                              %------------------                                                             %------------------
end                                                                                                                      %==================

end

%-----------------------------------------------------------------------------------------------------------------------

function [X,EM] = EMimputation(X,q0,opts)
% ===========================================================================================
% Function to impute missing values in a panel by means of a static factor
% model. Adapted from McCracken and Ng (2016).
% ------------------------------------------------------------------------
% SYNTAX: EMimputation(X,opts)
% where:
%      - X = matrix of data (TxN)
%      - q0 = number of factors (for testing)
%             if q0==99, select n. factors with Bai&Ng IC
%      - opts = structure containing:
%               q = number of factors to be estimated
%               maxiter = max. number of EM iterations (optional)
%               thresh = threshold to stop the algorithm (optional)
%               algo = algorithm for imputation ('SW' or 'BM')
%               T19 = end of pre-Covid period (for VAR)
% -------------------------------------------------------------------------
% OUTPUT: X, pc (structure)
% where:
%      - X = data with imputed values for NaNs (TxN)
%      - EM = structure containing:
%             F = matrix of estimated factors (Txq)
%             C = matrix of estimated loadings (Nxq)
%             chi = estimated common component (TxN)
%             R = estimated idosyncratic variances (TxN,diagonal)
%             d = q-largest eigenvalues (qx1)
% ------------------------------------------------------------------------------------------
% REFERENCES:
%
% [1] McCracken, M.W., and Ng, S. "FRED-MD: A monthly database for macroeconomic research."
%     Journal of Business & Economic Statistics, 34(4), 574-589. (2016)                [MN16]
% [2] McCracken, M.W., and Ng, S. "FRED-QD: A quarterly database for macroeconomic research."
%     No. w26872. National Bureau of Economic Research, 2020. (2020)                   [MN20]
% [3] Stock, J.H., and Watson, M.W. "Macroeconomic forecasting using diffusion indexes."
%     Journal of Business & Economic Statistics, 20(2), 147-162. (2002)                [SW02]
% [4] Bańbura, M., and Modugno, M. (2014). Maximum likelihood estimation of factor 
%     models on datasets with arbitrary pattern of missing data. 
%     Journal of applied econometrics, 29(1), 133-160.                                 [BM14]
% -------------------------------------------------------------------------------------------
% Author: Claudio Lissona (claudio.lissona2@unibo.it)
% Last update: 25/Aug/2023
% Version: MATLAB 2021b
% Required Toolboxes: /
% ===========================================================================================

if ~isfield(opts,'maxiter'); opts.maxiter = 500; fprintf('Maximum number of iterations set to %d',opts.maxiter); end     % default max-iterations: 500
if ~isfield(opts,'thresh'); opts.thresh = 10^-4; fprintf('Threshold for EM objective set to %10f \n',opts.thresh); end   % default threshold: 10e-4
if ~isfield(opts,'out'); opts.out = 0;  end                                                                              % default: no outlier adjustment

if opts.out~=0; [X,~,~] = remove_outliers(X); end                                                                        % remove outliers (if specified)
[T,N] = size(X);

indNaN = isnan(X);                                                                                                       % locate NaNs
mx = repmat(mean(X,'omitnan'),size(X,1),1);                                                                              % data mean, no NaNs
sx = repmat(std(X,'omitnan'),size(X,1),1);                                                                               % data std, no NaNs

Xnan = (X - mx)./sx;                                                                                                     % std. data w/ NaNs
X_est = Xnan; X_est(indNaN) = 0;                                                                                         % std. data w/out NaNs

X(indNaN) = mx(indNaN);                                                                                                  % impute missing values with sample mean

snan = sum(indNaN,2);                                                                                                    % for rows w/ more than 75% NaNs...
loc = find(snan>0.75*N);                                                                                                 % ...find rows...
for nn=1:numel(loc); Xnan(loc(nn),indNaN(loc(nn),:)) = 0; end                                                            % ...and set to unconditional mean (no imputation)

Xs = (X-mx)./sx;                                                                                                         % std. data w/out NaNs

if q0==99                                                                                                                % if n. of factors to be selected    
 if strcmp(opts.algo,'SW'); q = BaiNg(Xs,15,2); else; q = BaiNg(Xs,6,2); end                                             % set max to 15 for SW and to 6 for BM            
 else; q = q0;                                                                                                           % otherwise, set to q0          
end 

pc = princfact(Xs,q);                                                                                                    % PCA to extract factors
chi0 = pc.chi;                                                                                                           % store common component

ss = buildStateSpace(pc,opts);                                                                                           % build state-space from PC
EM = ss; EM.chi = pc.chi;                                                                                                % initialize EM output                   

KF = kalman(Xnan,ss);                                                                                                    % Kalman Filter/Smoother
FtT = KF.ZtT; PtT = KF.PtT; PtTm = KF.PtTm;                                                                              % extract factors and cond. var-cov
logl0 = KF.loglik;																																																																																														        	% initial log-likelihood 


err = 999; jj = 0;                                                                                                       % initialize EM 
%%%%%%%%%%%%%%%%%%%%
%%% EM ALGORITHM %%%
%%%%%%%%%%%%%%%%%%%%
while err>opts.thresh&&jj<opts.maxiter
    if mod(jj,10)==0
        fprintf('Running iteration %d: error is %10f. \n',jj,err);                                                       % updates every ten iterations
    end

    %%%%%%%%%%%%%%%%%%
    %%% IMPUTATION %%%
    %%%%%%%%%%%%%%%%%%                                                                                                   %=================================%
    if strcmp(opts.algo,'SW')                                                                                            %     1. STOCK & WATSON (2002)    %
       %---------------------                                                                                            %=================================%
       X(indNaN) = pc.chi(indNaN).*sx(indNaN) + mx(indNaN);                                                              % impute NaNs with common component
   
       mx = repmat(mean(X),size(X,1),1); sx = repmat(std(X),size(X,1),1); Xs = (X-mx)./sx;                               % standardize (new) panel
   
       if q0==99; q = BaiNg(Xs,15,2); else; q = q0; end                                                                  % select n. of factors
       
       pc = princfact(Xs,q); EM.chi = pc.chi;                                                                            % PCA to extract factors
   
       dchi = EM.chi-chi0; err = (dchi(:)'*dchi(:))/(chi0(:)'*chi0(:));                                                  % compute objective function ---> see [MW16]
   
       chi0 = EM.chi;                                                                                                    % update values for next iteration
       %---------------------                                                                                            %=================================%
    elseif strcmp(opts.algo,'BM')                                                                                        %   2. BANBURA & MODUGNO (2014)   %
       %---------------------                                                                                            %=================================%

       EM.Z00 = KF.ZtT(1,:)'; EM.P00 = KF.PtT(:,:,1);                                                                    % state initial values 

       FtFt1 = (FtT(2:opts.T19,:)'*FtT(1:opts.T19-1,:) + sum(PtTm(:,:,2:opts.T19),3));                                   % \sum_{t=1}^T F_{t}F_{t-1}'
       Ft1Ft1 = (FtT(1:opts.T19-1,:)'*FtT(1:opts.T19-1,:) + sum(PtT(:,:,1:opts.T19-1),3));                               % \sum_{t=1}^T F_{t-1}F_{t-1}'
       FtFt = (FtT(2:opts.T19,:)'*FtT(2:opts.T19,:) + sum(PtT(:,:,2:opts.T19),3));                                       % \sum_{t=1}^T F_{t}F_{t}'

       EM.C = (X_est(1:opts.T19,:)'*FtT(1:opts.T19,:))/FtFt;                                                             % estimate loadings

       EM.A = FtFt1/Ft1Ft1;                                                                                              % estimate state transition matrix
       EM.Q = (FtFt - EM.A*FtFt1')/T;                                                                                    % estimate state error var-cov matrix
       
       R0 = zeros(N,N);                                                                                                  % initialize idiosyncratic var-cov matrix
       for tt = 1:opts.T19
           Wt = diag(~indNaN(tt,:));                                                                                     % selection matrix: non-missing values at time t
           R0 = R0 + (X_est(tt,:)' - Wt*EM.C*FtT(tt,:)')*(X_est(tt,:)' - Wt*EM.C*FtT(tt,:)')'...
                 + Wt*EM.C*PtT(:,:,tt)*EM.C'*Wt + (eye(N)-Wt)*EM.R*(eye(N)-Wt);                                          % [BM14] eq. 12
       end
       EM.R = diag(diag(R0)/T);                                                                                          % ver-cov matrix idiosyncratic component 

       EM.chi = FtT*EM.C';                                                                                               % common component                                                         
        
       KF = kalman(Xnan,EM);                                                                                             % Kalman Filter/Smoother
       FtT = KF.ZtT; PtT = KF.PtT; PtTm = KF.PtTm;                                                                       % update states

       logl = KF.loglik;																																																																																															  % new log-likelihood     	
       err = (abs(logl - logl0)/((abs(logl) + abs(logl0))/2));																																																		         % --------------
    			logl0 = logl;																																																																																																	    % update previous likelihood
       %---------------------
    end
    jj = jj+1;
    if err<opts.thresh&& jj<opts.maxiter; fprintf('EM converged after %d iterations \n',jj); end                         % show convergence (if any)

X(indNaN) = EM.chi(indNaN).*sx(indNaN) + mx(indNaN);                                                                     % impute NaNs with common component

end

end


function ss = buildStateSpace(pc,opts)

p = 1;                                                                                                                   % fix VAR lag to 1
q = size(pc.F,2);                                                                                                        % get number of factors                            

res = estVAR(pc.F(1:opts.T19,:),p);                                                                                      % estimate VAR       

ss.A = res.J; ss.Q = zeros(q,q); ss.Q(1:q,1:q) = res.Q;                                                                  % parameters: state eq.
ss.C = pc.C; ss.R = diag(diag(pc.R));                                                                                    % parameters: observation eq.
ss.Z00 = []; for pp=1:p; ss.Z00=cat(1,ss.Z00,pc.F(p+1-pp,:)'); end                                                       % initial states: factors
ss.P00 = reshape((eye((q)^2)-kron(ss.A,ss.A))\ss.Q(:), q, q);                                                            % initial var-cov states: factors 

end