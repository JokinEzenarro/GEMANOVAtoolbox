function model = gemanova(X,F,Fix,scl,cross,show)

%GEMANOVA - GEneralized Multiplicative ANOVA
% Fits a GEMANOVA model to an N-way array of responses 
% X with F effects. 
%
% INPUTS
% The third input, Fix defines the nature 
% of the effects. Fix is a binary matrix of size 
% N x F and the f'th column defines the f'th effect. 
% If the column contain zeros only, the effect is
% a multiplicative effect of all modes.
%
% Ex.: For a three-way array of size I x J x K, 
%
% Fix(:,f) = [0 0 0]' corresponds to a_i*b_j*c_k
% Fix(:,f) = [1 0 0]' corresponds to b_j*c_k
% Fix(:,f) = [1 0 1]' corresponds to b_j
% Fix(:,f) = [1 1 1]' corresponds to m (an overall constant)
% 
% Thus a two-effect model: x_ijk = a_i*b_j*c_k + c_k
% can be specified as
%
% model = gemanova(X,2,[0 0 0;1 1 0]',scl);
%
% scl is an optional input of a cell array of vectors for 
% plotting effects against (set scl=[] for none)
%
% A fifth optional input, cross, performs cross-validation
% if set to one. 
% 
% OUPUTS
% The output is a structured array model which holds
%     xname: name of the original workspace input variable
%      name: type of model, always 'GEMANOVA'
%      date: model creation date stamp
%      time: model creation time stamp
%      size: size of the original input array
%  neffects: number of effects estimated
%       ssq: residual sum of squares
%   effects: 1 by order cell array of the loadings in each dimension
%       res: 1 by order cell array squared residuals summed over each dimension
%       scl: 1 by order cell array with scales for plotting loads 
%
% If cross-validation is performed, the additional output
% is a structure Cross which holds the reproduced 
% array obtained by leave-one-element-out cross-validation 
% (as the elements are assumed to be independent in ANOVA) 
% as well as the RMSECV and the correlation between the 
% responses and the reproduced array. Note that MANOVA can be
% performed directly with gemanova but the cross-validation
% assumes independent elements.
%
% The algorithm handles missing elements if these are 
% set to NaN of Inf. 
%
% Please refer to this m-file through
%
%      Rasmus Bro & Marianne Jakobsen, Exploring complex 
%      interactions in designed data using GEMANOVA. Color 
%      changes in fresh beef during storage, Journal of 
%      Chemometrics, 2001, Submitted.
%
% I/O: model = gemanova(X,2,[0 0 0;1 1 0]',scl,Cross);

%rb February, 2001

id = find(isinf(X));
X(id) = NaN; % Because the following algorithms only handle NaN

if (nargin < 4 | ~strcmp(class(scl),'cell'))
    scl = cell(1,length(size(X)));
end
if nargin < 5
    cross = 0;
end
if nargin < 6
    show = 1;
end
if any(sum(Fix)==length(size(X)));
    % Means that there is a constant across all modes
    % Then an extra mode of dimension one has to be added
    % in the end, the mode will be removed again
    x = zeros([1 size(X)]);
    x(1,:) = X(:)';
    X = x;
    clear x;
    Fix = [ones(1,F);Fix];
    i = find(sum(Fix)==length(size(X)));
    Fix(1,i)=0;
    ExtraMode = 1;
else
    ExtraMode = 0;
end

if cross==1
    % Find overall model
    disp(' ')
    disp(' CROSS-VALIDATING THE SOUGHT MODEL')
    disp('  Fitting the global model ...')
    model = gemanova(X,F,Fix,scl);
    disp(' ')
    disp('  Cross-validating ...')
    Cross = gemcross(X,F,Fix,scl,model.nparam);
    model.Cross = Cross;  
    
else % Fit model 
    if show
        disp(' ')
        disp(' Gradually adjusting to correct model')
    end
    
    model = parafacweight(X,F,Fix,show);
    if show
        disp(' Refining to exact solution')
    end
    [m,fit]=parafacones(X,F,Fix,5000,model,show);
    % fix added mode if necessary
    if ExtraMode
        Fix = Fix(2:end,:);
        clear effect
        for i = 1:length(m)-1
            effect{i} = m{i+1};
        end
        effect{end}=effect{end}*diag(m{1});
        m = effect;
        ExtraMode = 0;
        X = squeeze(X);
        if show
            disp(' Offsets constant across all modes are incorporated into the last mode parameters')
        end
        
    end
    
    
    % OUTPUT MODEL
    model = struct('xname',inputname(1),'name','GEMANOVA','date',date,'time',clock,...
        'size',size(X),'neffects',F);
    model.effects = m;
    model.ssq = fit;
    model.scl = scl;
    model.res = X-outerm(model.effects);
    nparam = 0;
    for i = 1:size(Fix,1)
        nparam = nparam + sum(Fix(i,:)==0)*size(X,i);
    end
    model.nparam = nparam;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% AUXILIARY FILES%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cross = gemcross(X,F,Fix,scl,nparam);

% Do cross-validation
Idx = find(~isnan(X));
Pred = NaN*ones(size(X));
disp([' ',num2str(length(Idx)),' segments']);
fprintf([' '])
FIT = [];
for i = 1:length(Idx)
    x = X;
    x(Idx(i)) = NaN;
    crossmodel = gemanova(x,F,Fix,scl,0,0);
    crossmod = outerm(crossmodel.effects);
    Pred(Idx(i))=crossmod(Idx(i));
    e = X(find(~isnan(x)))-crossmod(find(~isnan(x)));
    fit = sum(e(:).^2);
    FIT = [FIT;fit]; % Save this to check that there are no strange fit values.
    if i~=length(Idx)
        fprintf([num2str(i),','])
    else
        fprintf(num2str(i))
    end
%    disp([' Crossval segment ',num2str(i),' of ',num2str(length(Idx)),' - ssq: ',num2str(sum(e(:).^2))]);
end
Cross.Predict = Pred;
Cross.RMSECV = std(X(Idx)-Pred(Idx));
e = (X(Idx)-Pred(Idx)).^2;
Cross.RMSECVcorrect = sqrt( sum(e)/(length(Idx)-nparam));
r = corrcoef([vec(X(Idx)) vec(Pred(Idx))]);
Cross.corr = r(2,1);
Cross.IndividualSsq = FIT;




function model = parafacweight(X,F,Fix,show)

% Parafac imposing fixed elements with ones gradually through weighting
% The end result will not have exact ones but will be a good initial 
% starting point for an exact algoritm


crit = 1e-10;
showfit = pi; % Don't show anything
showfitinit = pi;
numsplitgradient = 6;
maxit = 500;
if nargin<4
    show = 1;
end
order = length(size(X));
it = 0;
if size(Fix,1)~=order & size(Fix,2)~=F
    error(' Fix is not correctly given')
end
DimX = size(X);

% Initialize with parafac
fit = sum(X(find(~isnan(X))).^2);
fitold = fit*2;
m = parafac(X,F,[],[],[],[1 0],0);
for rep = 1:3 % See if randomly started is better
    x0 = m.loads;
    for i = 1:length(x0)
        x0{i} = rand(size(x0{i}));
    end
    m2 = parafac(X,F,[],[],x0,[1 0],0);
    if m2.ssq(2)<m.ssq(2)
        m = m2;
    end
end        
loads = m.loads;

% Find gradually increasing weights
ssqX = fit;
Weights = linspace(0,ssqX/300,numsplitgradient);

for w = 1:numsplitgradient
    fitold = fit*10;
    it = 0;
    while abs(fit-fitold)/fitold>crit & it < maxit
        
        it = it+1;
        fitold = fit;
        
        % Adjust scales on different loads to avoid numerical problems
        scales = ones(F,1);
        for fac = 1:F
            for ord = 1:order
                if ~Fix(ord,fac)
                    scales(fac) = scales(fac)*norm(loads{ord}(:,fac));
                end
            end
        end
        for fac = 1:F
            NumbNonFix = length(find(~Fix(:,fac)));
            scales(fac) = scales(fac)^(1/NumbNonFix);
            for ord = 1:order
                if ~Fix(ord,fac)
                    loads{ord}(:,fac) = scales(fac)*loads{ord}(:,fac)/norm(loads{ord}(:,fac));
                end
            end
        end
        
        for fac = 1:F
            for ord = 1:order
                % find update for ord'th mode of fac'th factor
                z = outerm(loads,ord,1);
                permX = permute(X,[ord 1:ord-1 ord+1:order]);
                permXunf = reshape(permX,DimX(ord),prod(DimX([1:ord-1 ord+1:end])));
                % handle missing
                for p = 1:DimX(ord)
                    id = find(~isnan(permXunf(p,:)));
                    xx= permXunf(p,id)';
                    zz = z(id,:);
                    if Weights(w)&any(Fix(ord,:))
                        xx = [xx;ones(length(find(Fix(ord,:))),1)*Weights(w)];
                        ffix = eye(F);
                        ffix(find(~Fix(ord,:)),:)=[];
                        zz  = [zz;ffix*Weights(w)];
                    end
                    ll = xx'*pinv(zz)';
                    if sum(abs(ll))<100*eps; % almost zero
                        loads{ord}(p,:) = .9*loads{ord}(p,:)+.1*ll;
                    else
                        loads{ord}(p,:) = ll;
                    end
                end
            end
            
        end
        
        
        fit = pffit(X,loads);
        if show & rem(it,showfit)==0 & it > 1
            disp([' Fit of GEMANOVA model = ',num2str(fit),' after ',num2str(it),' it. in final round & Weight = ',num2str(Weights(w)),' (',num2str(w),'/',num2str(numsplitgradient),') '])
        end
    end
end

model = loads;
if show
    disp([' Final fit of GEMANOVA model = ',num2str(fit),' after ',num2str(it),' it. & Weight = ',num2str(Weights(w)),' (',num2str(w),'/',num2str(numsplitgradient),')'])
end
%%%%%%%%%%%%%%%%%%% END OF PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ssq = pffit(X,loads);
M = outerm(loads);
E = X-M;
ssq = sum(E(find(~isnan(E))).^2);


function product = khatri(loads,factor);
product = [];
for fac = 1:size(loads{1},2)
    if fac~=factor
        Z = kron(loads{end}(:,fac),loads{end-1}(:,fac));
        for j = length(loads)-2:-1:1
            Z = kron(Z,loads{j}(:,fac));
        end
        product = [product Z];
    end
end


function [m,fit]=parafacones(X,F,Fix,maxit,loads,show)

% Fit exact ls parafac with exact ones

crit = 1e-11;
showfit = pi;
showfitinit = pi;
if nargin < 4
    maxit = 1000;
end
if nargin<6
    show = 1;
end
order = length(size(X));
it = 0;
if size(Fix,1)~=order & size(Fix,2)~=F
    error(' Fix is not correctly given')
end
DimX = size(X);


fit = sum(X(find(~isnan(X))).^2);
fitold = fit*2;
while abs(fit-fitold)/fitold>crit & it < maxit
    
    it = it+1;
    fitold = fit;
    
    if it >5 & rem(it,10)==0 %& 5==3% Do line-search
        options = optimset;
        options = optimset(options,'Display','off');
        [alpha,fithere,exitf] = fminbnd('pffitalpha',1,100,options,X,loads,LoadingsOld);
        if fithere<fit
            for ord = 1:order
                loads{ord} = loads{ord}+alpha*(loads{ord}-LoadingsOld{ord});
            end
            fit = fithere;
        end
        if show & rem(it,showfit)==0
            disp([' Fit after line-search = ',num2str(fit),' after ',num2str(it),' iterations'])
        end
    end
    
    LoadingsOld = loads;
    for fac = 1:F
        
        % Calculate residual of X for model with all but fac'th factor
        if F>1
            product = khatri(loads,fac);
        else
            product = zeros(size(X(:)));
        end
        if F>2
            sumprod = sum(product');
        else
            sumprod = product;
        end
        resX = X - reshape(sumprod,size(X));
        
        % calculate old loads for fac
        for ord = 1:order;
            oldloads{ord} = loads{ord}(:,fac);
        end
        
        
        for ord = 1:order
            % find update for ord'th mode of fac'th factor
            z = outerm(oldloads,ord,1);
            if ~Fix(ord,fac)
                permresX = permute(resX,[ord 1:ord-1 ord+1:order]);
                permresXunf = reshape(permresX,DimX(ord),prod(DimX([1:ord-1 ord+1:end])));
                % handle missing
                for p = 1:DimX(ord)
                    id = find(~isnan(permresXunf(p,:)));
                    oldloads{ord}(p) = permresXunf(p,id)*pinv(z(id))';
                end
            else
                oldloads{ord}=ones(DimX(ord),1);
            end
        end
        
        % update loads
        for ord = 1:order
            loads{ord}(:,fac) = oldloads{ord};
         end
         
    end
    
    
    fit = pffit(X,loads);
    if show & rem(it,showfit)==0 & it > 1
        disp([' Fit of GEMANOVA model = ',num2str(fit),' after ',num2str(it),' iterations'])
     end
end

% Scale all but last (non-fixed) mode to unit length
for f = 1:F
    for o = 1:order-1
        if ~Fix(o,f) % Not fixed to ones, hence should be normalized
           if any(~Fix(o+1:end,f)) % Only normalize if there is a 'later' unfixed mode
              oo = o+1:order;
                fix = find(~Fix(oo,f));
                fix = oo(fix(end));
                loads{fix}(:,f) = loads{fix}(:,f)*norm(loads{o}(:,f));
                loads{o}(:,f) = loads{o}(:,f)/norm(loads{o}(:,f));
            end
        end
    end
end

m = loads;
if show
    disp([' Final fit of GEMANOVA model = ',num2str(fit),' after ',num2str(it),' iterations'])
 end
 
 
 function vX = vec(X);
 
 vX = X(:);
 
 
 
 
 
 
 
 
 
 
 function model = parafac(x,nocomp,scl,tol,x0,options,plots,constraints,weights);

%PARAFAC Parallel factor analysis for n-way arrays
%  PARAFAC will decompose an array of order n (where n >= 3)
%  into the summation over the outer product of n vectors.
%  Missing values must be NaN or Inf.
%  
%  INPUTS
%  x           : the multi-way array to be decomposed
%  nocomp      : the number of components to estimate, 
%
%  OPTIONAL INPUTS
%  scl         : optional inputs of a cell array of vectors for plotting
%                loads against (set scl=[] for none, if elements of scl are
%                empty or of incorrect length they are reset),
%  tol         : convergence tolerance (a 1 by 3 vector consisting of 
%                the relative change in fit, absolute change in fit 
%                and maximum iterations {default tol = [1e-6 1e-6 10000]}), 
%  x0          : initial estimate of the factors, default is none = [].
%                x0 is a cell array of the same form as the fitted loadings
%                output in model.loads. Thus, x0{i} holds the loading matrix
%                in the i'th mode.
%                   If x0 is correctly provided, it will be used for initializing
%                the algorithm. Note, that if constraints(i) is set to -1, then
%                the loadings of that mode will not be updated, hence the initial
%                stay fixed. This way, one may use the algorithm for fitting an 
%                existing model to new data, by inputting the old model (model.loads) 
%                and fixing all except the sample mode loadings.
%                   It is also possible to input a complete PARAFAC model (a 
%                structure). This is useful when a prior model also includes
%                offsets, so that these are also used as initial values.
%  options     : A vector defining algorithmic settings:
%                options(1) = a flag which turns off the line search options when 0
%                options(2) = turns off screen output when set to zero
%                Default, options = [] = [1 1];
%  plots       : a flag which turns off the plotting of the loads and 
%                residuals when set to zero. 
%  constraints : Vector defining constraints. constraints has length equal to the number
%                of modes (for a three-way array, constraints is a 3-vector). 
%                constraints(1) =  0 => First mode loadings unconstrained
%                constraints(1) =  1 => First mode loadings nonnegativity-constrained
%                constraints(1) =  2 => First mode loadings unimodality-unconstrained (and nonneg)
%                constraints(1) =  3 => First mode loadings orthogonality-unconstrained
%                constraints(1) = -1 => First mode loadings fixed (requires initial values given)
%                constraints(2) =  ...  Defines constraints on second mode loadings
%                Example.: constraints = [0 1 2] means A unconstrained, B nonneg, & C unimodal
%  weights     : Enables weighted least squares fitting. Type <parafac('weights')> to get help on 
%                this feature
%
%  The output is a structured array (model) containing the 
%  PARAFAC model elements:
%
%     xname: name of the original workspace input variable
%      name: type of model, always 'PARAFAC'
%      date: model creation date stamp
%      time: model creation time stamp
%      size: size of the original input array
%     ncomp: number of components estimated
%       tol: tolerance vector used during model creation
%     final: final tolerances at termination
%       ssq: total and residual sum of squares
%     loads: 1 by order cell array of the loadings in each dimension
%       res: 1 by order cell array squared residuals summed over each dimension
%       scl: 1 by order cell array with scales for plotting loads 
%
%  This routine uses alternating least squares (ALS) in combination with
%  a line search every fifth iteration. For 3-way data, the intial estimate
%  of the loadings is obtained from the tri-linear decomposition (TLD).
%
%I/O: model = parafac(x,nocomp,scl,tol,x0,options,plots,constraints,weights,iteroff);
%
%See also:  GRAM, MPCA, MWFIT, OUTER, OUTERM, TLD, UNFOLDM, XPLDST

if nargin == 0
	disp(' ')
	disp(' PARAFAC Parallel factor analysis for n-way arrays - short help')
	disp(' ')
	disp(' I/O: model = parafac(x,nocomp,scl,tol,x0,options,plots,constraints,weights,iteroff);')
	disp(' ')
	disp(' Type <<help parafac>> for extended help')
	disp(' ')
	return
end

if nargin ==1
	if strcmp(lower(x),'weights')
		disp(' USING A WEIGHTED LOSS FUNCTION')
		disp(' ')
		disp(' Through the use of the input ''weights'' it is possible to fit a PARAFAC')
		disp(' model in a weighted least squares sense')
		disp(' ')
		disp(' The input is an array of the same size as the input data holding individual')
		disp(' weights for changing the loss function from a least squares to weighted least')
		disp(' squares one. Instead of minimizing the frobenius norm ||x-M|| where M is the PARAFAC')
		disp(' model, the norm ||(x-M).*weights|| is minimized. The algorithm used for')
		disp(' weighted regression is a majorization step according to Kiers, Psychometrika,')
		disp(' 1997, 62, 251-266, which has the advantage of being computationally cheap and ')
		disp(' also handles situations where ordinary weighted least squares regression does')
		disp(' not work.')
	end
	return
end



% INITIALIZATION OF ALGORITHM
Show = 100; % How often are fit values shown
xsize = size(x);
order = ndims(x);
if (nargin < 3 | ~strcmp(class(scl),'cell'))
	scl = cell(1,order);
end
for ii=1:order   %added 3/1/00 nbg
	if isempty(scl{ii})|(length(scl{ii})~=xsize(ii))
		scl{ii} = 1:xsize(ii);
	end
end
standardtol = [1e-6 1e-6 10000];
if (nargin < 4 | length(tol)==0)
	tol = standardtol;
else % Changed error rb 11/03/00
	tol(find(tol==0)) = standardtol(find(tol==0));
end
if nargin < 6
	options = [1 1];
end
if length(options)<2
	options = [options ones(1,2-length(options))];
end
ls = options(1); %For line-searh
DumpToScreen = options(2);
if nargin < 7|length(plots)==0
	plots = 1;
end
if nargin < 8
	constraints = zeros(1,order);
end
if length(constraints)~=order
	constraints = zeros(1,order);
end

if nargin < 9 
	DoWeight = 0;
else
	if length(size(weights))~=order 
		DoWeight = 0;
	elseif all(size(weights)==xsize)
		DoWeight = 1;
		WMax     = max(abs(weights(:)));
		W2       = weights.*weights;
	else
		DoWeight = 0;
	end
end


%Iterative preproc
iteroffsets = 0;
if exist('iteroff')==1
	if iscell(iteroff)
		iteroffsets = 1;
	elseif isnumeric(iteroff) % Then convert to cell
		hh = iteroff;
		clear iteroff
		iteroff{1} = hh;
		iteroffsets = 1;
	elseif ischar(iteroff)
		if strcmp(lower(iteroff(1:2)),'he')
			clear iteroff
			button = questdlg('Do you want help to incorporate offsets?','Offsets','Yes','No','Cancel','Yes');
			if strcmp(button,'Cancel')
				return
			elseif strcmp(button,'No')
				iteroffsets = 0;
			else
				Continue = 1;
				while Continue
					if Continue==1
						questdlg('For each type of offset, you input the number of the modes across which the offsets are constant. E.g. for a three-way array where you want ordinary offsets that are constant across the first mode (~average of each column); then you input the offset as [1], meaning the the offset is constant across the first mode. In this case the offsets will therefore be a matrix of dimension J x K (J second mode dimension, K third mode dimension). If you want an offset that is constant across all three modes, you input [1 2 3], implying that one offset is estimated','Offsets - how to do it','OK','OK');
					end
					prompt={['Enter the definition of offset # ',num2str(Continue)]};
					def={'[1]'};
					dlgTitle='Offset definition';
					lineNo=1;
					% Ask for offset
					Def=inputdlg(prompt,dlgTitle,lineNo,def);
					% If cancel is hit instead of OK
					if isempty(Def)
						hh=questdlg('Offsets across the first mode will be assumed','Drop offsets?','OK','Stop PARAFAC','OK');
						if strcmp(hh,'Stop PARAFAC')
							error(' PARAFAC stopped')
						else
							Def{1} = '[1]';
						end
					end
					% Verify the offset and show how it's implemented
					hh = eval(Def{1});
					txt = ['[',num2str(hh(1))];
					for i = 2:length(hh);
						txt = [txt,' ',num2str(hh(i))];
					end
					txt = [txt,']'];
					prompt={['Offset # ',num2str(Continue),' is defined as: iteroff{',num2str(Continue),'} = ',txt,';']};
					dlgTitle='Learn how to use the offset input <iteroff>';
					% Ask for offset
					questdlg(prompt,dlgTitle,'OK','Cancel','OK');
					iteroff{Continue} = eval(Def{1});
					Continue = Continue + 1;
					answer=questdlg('Do you want to incorporate other offsets?','More?','Yes','No','Yes');
					if strcmp(answer,'No')
						Continue = 0;
					end
				end
				iteroffsets = 1;
			end
		end
	end
end
if iteroffsets
	XTrilin = zeros(size(x));
	%   options(1) = 0; % No line-search can be performed as iterative offsets are not currently supported in the line-search
end

FeasibilityProblems = 0; % for indicating if the algorithm hasn't yet reached a feasible solution due some constraints
%Define in which mode the scale should go - default is last unless it's fixed')
if all(constraints==-1)
	error(' All modes are fixed, hence no fitting problem')
else
	if constraints(order)==-1 % Since last mode is fixed, scale of the components cannot go there. Instead use the first non-fixed mode
		ScaleMode = min(find(constraints~=-1));
	else
		ScaleMode = order;
	end
end

% CHECK FOR MISSING
if any(isinf(x(:)))|any(isnan(x(:)))
	MissId    = sparse(find(isinf(x)|isnan(x)));
	Missing   = 1;
else
	Missing = 0;
end

% INITIALIZE LOADINGS (if not already given)
if (nargin < 5 | ~(strcmp(class(x0),'cell') | strcmp(class(x0),'struct')))  % Old loadings not given
	x0 = cell(1,order);
	if all(xsize>nocomp) & ~Missing% Use atld
		x0=atld(x,nocomp,0);
		InitString = ' Using fast approximation for initialization';
	elseif order == 3&~Missing&~all(xsize<nocomp)
		% Initialize with TLD estimates
		m = tld(x,nocomp,0,0);
		x0 = m.loads;
		InitString = ' Using direct trilinear decomposition for initialization';
  else
		for j = 1:order
			x0{1,j} = rand(xsize(j),nocomp);
		end
		InitString = ' Using random values for initialization';
	end
else % Old loadings given
	if strcmp(class(x0),'cell')
		x0 = x0;
		InitString = ' Using old values for initialization';
	elseif strcmp(class(x0),'struct')
		if ~strcmp(x0.name,'PARAFAC')
			error(' Input x0 is model that is not a PARAFAC model (name in structure should be PARAFAC)')
		else
			% Use prior model given in x0 to extract loadings and possibly offsets
			if isfield(x0,'Offsets')  % If offsets are given in prior model, then include these
				iteroffsets = 1;
				iteroff = x0.Offsets.Def;
				for i=1:length(iteroff);
					Offset{i}.OffsetsReduced = x0.Offsets.Parameters{1};
					Offset{i}.Offsets = x0.Offsets.Model{1};
				end
				x0 = x0.loads;
				XTrilin = outerm(x0);%zeros(size(x));
				InitString = ' Using old values for initialization (including offsets)';
			else
				x0 = x0.loads;
				InitString = ' Using old values for initialization';
			end
		end
	elseif ~strcmp(class(x0),'cell')
		error('Initial estimate x0 not a cell array')
	end
end

% CHECK FOR INFEASIBLE INITIAL SOLUTIONS
if any(constraints)==1|any(constraints==2) % Nonnegativity required => make sure initial values are positive
	for j=1:order
		if constraints(j)==1|constraints(j)==2
			if any(x0{j}(:)<0) % Negative values occur
				x0{j} = sign(sum(sum(x0{j})))*x0{j}; %'Turn' it around so mainly positive
				i=find(x0{j}<0);
				x0{j}(i)=-x0{j}(i);
			end
		end
	end
end
if any(constraints==3) % Orth required
	for j=1:order
		if constraints(j)==3
			if size(x,j)<nocomp
				error([' Orthogonality cannot be applied in a mode (mode ',num2str(j),') where the dimension is less than the number of factors'])
			else
				x0{j} = orth(x0{j});
			end
		end
	end
end

if Missing % exchange missing elements with model estimates
	xest = outerm(x0);   % Old version;xest = zeros(prod(xsize),nocomp);   for j = 1:nocomp      xvect = x0{1}(:,j);      for ii = 2:order         xvect = xvect*x0{ii}(:,j)';         xvect = xvect(:);      end      xest(:,j) = xvect;   end,   xest = sum(xest,2);   xest = reshape(xest,xsize);
	XTrilin = xest; % Multilinear part to subtract from data before calculating offsets
	if iteroffsets & exist('Offset') % it doesn't exist unless prior model is input
		for i=1:length(iteroff) % estimate each set of offsets one at a time
			xest = xest + Offset{i}.Offsets;
		end
	end
	x(MissId)=xest(MissId);
end
% Calculate total sum of squares in data
if DoWeight
	xsq = (x.*weights).^2;
else
	xsq = x.^2;
end
if Missing
	xsq(MissId)=0;  
end
tssq = sum(xsq(:));

% Initialize the unfolded matrices
xuf = cell(1,order);
xuflo = cell(1,order);
xufsize = zeros(1,order);
for i = 1:order
	xuf{i} = unfoldmw(x,i)';
	xufsize(i) = prod(xsize)/xsize(i);
	% Old init of loads;xuflo{i} = zeros(prod(xsize)/xsize(i),nocomp);
	for j = 1:nocomp
		% Old version if i == 1,         mwvect = x0{2}(:,j);,         for k = 3:order              mwvect = mwvect*x0{k}(:,j)';            mwvect = mwvect(:);         end      else         mwvect = x0{1}(:,j);         for k = 2:order            if k ~= i               mwvect = mwvect*x0{k}(:,j)';               mwvect = mwvect(:);            end         end      end      xuflo{i}(:,j) = mwvect;   end
		xuflo{j}=outerm(x0,j,1);
	end
end

% Initialize other variables needed in the ALS
oldx0 = x0;
searchdir = x0;
iter = 0;
flag = 0;
oldests = zeros(prod(xsize)*nocomp,1);
abschange = 0;
relchange = 0;


% Show algorithmic settings etc.
if DumpToScreen
	disp(' ')
	disp(' Fitting PARAFAC ...') 
	txt=[];
	for i=1:order-1
		txt=[txt num2str(xsize(i)) ' x '];
	end
	txt=[txt num2str(xsize(order))];
	disp([' Input: ',num2str(order),'-way ',txt, ' array'])
	
	disp([' A ',num2str(nocomp),'-component model will be fitted'])
	for i=1:order
		if constraints(i)==0
			disp([' No constraints on mode ',num2str(i)])
		elseif constraints(i)==1
			disp([' Nonnegativity on mode ',num2str(i)])
		elseif constraints(i)==2
			disp([' Unimodality on mode ',num2str(i)])
		elseif constraints(i)==3
			disp([' Orthogonality on mode ',num2str(i)])
		elseif constraints(i)==-1
			disp([' Fixed loadings in mode ',num2str(i)])
		end
	end
	
	disp([' Convergence criteria:']) 
	disp([' Relative change in fit : ',num2str(tol(1))]) 
	disp([' Absolute change in fit : ',num2str(tol(2))]) 
	disp([' Maximum iterations     : ',num2str(tol(3))]) 
	
	if Missing
		disp([' ', num2str(100*(length(MissId)/prod(xsize))),'% missing values']);
	else
		disp(' No missing values')
	end
	
	if DoWeight
		disp(' Weighted optimization will be performed using input weights')
	end
	
	%Iterative preproc
	if iteroffsets
		disp(' Iteratively calculated offsets will be used')
		disp(' Line-search has been turned off')
	else
		%disp(' No iterative offsets')
	end
	
	disp(InitString)
end

% Start the ALS
while flag == 0;
	iter = iter+1;
	% Loop over each of the order to estimate
	
	if DoWeight & iter > 1% If Weighted regression is to be used, do majorization to make a transformed data array to be fitted in a least squares sense
		out = reshape(xest,xsize) + (WMax^(-2)*W2).*(x - reshape(xest,xsize));
		for i = 1:order
			xuf{i} = unfoldmw(out,i)';
		end
		clear out
	end
	
	%Iterative preproc
	if iteroffsets
		if iter>1
			OldOffset = Offset;
		end
		if DoWeight & Iter > 1 % Then the data to fit is not the raw data in x but the modified data according to above in xuf
			xprep = reshape(xuf{1}',xsize) - Xtrilin;
		else
			xprep = x - XTrilin;
		end
		for i=1:length(iteroff) % estimate each set of offsets one at a time
			%[OffsetsFromNstatReduced,OffsetsFromNstat] = nstat(xprep,'mean',iteroff{i},1);
			[OffsetsFromNstatReduced,OffsetsFromNstat] = nstat(xprep,'mean',iteroff{i});
			Offset{i}.OffsetsReduced = OffsetsFromNstatReduced;
			Offset{i}.Offsets = OffsetsFromNstat;
			xprep = xprep - OffsetsFromNstat;
		end
		% Replace unfolded data (xuf), so that the multilinear part is fitted to the residuals of the data subtracted the offsets
		if DoWeight & Iter > 1 % Then the data to fit is not the raw data in x but the modified data according to above in xuf
			xprep = reshape(xuf{1}',xsize);
		else
			xprep = x;
		end
		for i=1:length(iteroff) % estimate each set of offsets one at a time
			xprep = xprep - Offset{i}.Offsets;
		end
		for i = 1:order
			xuf{i} = unfoldmw(xprep,i)';
		end
	end
	
	for i = 1:order
		% Multiply the loads of all the orders together
		% except for the order to be estimated
		xuflo{i}=outerm(x0,i,1);
		
		% Regress the actual data on the estimate to get new loads in order i
		
		% Unconstrained
		if constraints(i) == 0
			ordiest = xuflo{i}\xuf{i};
			
			% nonnegativity   
		elseif constraints(i) == 1 
			ordiest = zeros(nocomp,xsize(i));
			xufloT = xuflo{i}'*xuflo{i}; % Calculate xproduct before loop
			for k = 1:xsize(i)
				ordiest(:,k) = CrossProdFastNnls(xufloT,xuflo{i}'*xuf{i}(:,k));
			end
			if any(sum(ordiest,2)==0);
				FeasibilityProblems=1;
				ordiest = .1*ordiest+.9*x0{i}';
			else
				FeasibilityProblems=0;
			end
			
			%Unimodality
		elseif constraints(i) == 2 
			ordiest=unimodal(xuflo{i},xuf{i},x0{i})';
			if any(sum(ordiest,2)==0);
				FeasibilityProblems=1;
				ordiest = .1*ordiest+.9*x0{i}';
			else
				FeasibilityProblems=0;
			end
			
			%Orthogonality
		elseif constraints(i) == 3
			% if this is the mode holding the scales modify so that loadings are not forced to be length one
			if i == ScaleMode
				Z = [];
				for fac = 1: nocomp
					Z = [Z kron(x0{i}(:,fac)/norm(x0{i}(:,fac)),xuflo{i}(:,fac))];
				end
				Scales = pinv(Z'*Z)*(Z'*xuf{i}(:));
			else
				Scales = ones(1,nocomp);
			end
			ZtX = (xuflo{i}*diag(Scales))'*xuf{i};
			ordiest=((ZtX*ZtX')^(-.5))*ZtX;
			if i == ScaleMode 
				ordiest = diag(Scales)*ordiest;
			end
			
			%Fixed
			
		elseif constraints(i) ~= -1 
			error([' The input constraints has not been correctly defined. The value ',num2str(constraints(i)),' is not possible'])
			
		elseif constraints(i) == -1 
			ordiest = x0{i}';
		end
		
		% Normalize the estimates (except the last (or other defined) order) and store them in the cell
		if i ~= ScaleMode & constraints(i)~=-1 % thus mode fixed
			Scal = 1./sqrt(sum(ordiest.^2,2));
			x0{i} = ordiest'*diag(Scal); % normalization leads to wrong model but that's corrected in the next update of the next mode, and for the last mode no normalization is performed, so that's ok, unless last mode is fixed.
			x0{ScaleMode} = x0{ScaleMode}*diag(Scal.^(-1));%ii = i+1;
		else
			x0{i} = ordiest';%ii = 1;
		end
		
	end
	% Calculate the estimate of the input array based on current loads
	xest = outerm(x0);
	% old version; xest = zeros(prod(xsize),nocomp);for j = 1:nocomp,      xvect = x0{1}(:,j);,      for ii = 2:order,         xvect = xvect*x0{ii}(:,j)';         xvect = xvect(:);      end,      xest(:,j) = xvect;,   end,   xest = sum(xest,2);   xest = reshape(xest,xsize);
	%Iterative preproc
	if iteroffsets
		XTrilin = xest; % Multilinear part to subtract from data before calculating offsets
		for ii=1:length(iteroff) % Add offsets one at a time
			xest = xest + Offset{ii}.Offsets;
		end      
	end
	xsq = xest;
	% Exchange missing with model estimates
	if Missing 
		x(MissId)=xest(MissId);
		for ii = 1:order
			xuf{ii} = unfoldmw(x,ii)';
		end
	end
	% Check to see if the fit has changed significantly
	if DoWeight
		xsq = ((x-xsq).*weights).^2;
	else
		xsq = (x-xsq).^2;
	end
	if Missing
		xsq(MissId)=0;  
	end
	ssq = sum(xsq(:));
	
	
	%disp(sprintf('On iteration %g ALS fit = %g',iter,ssq));
	if iter > 1 &~FeasibilityProblems
		abschange = abs(oldssq-ssq);
		relchange = abschange/ssq;
		if relchange < tol(1)
			flag = 1;
			if DumpToScreen
				disp(' '),disp('    Iteration    Rel. Change         Abs. Change         sum-sq residuals'),disp(' ')
				fprintf(' %9.0f       %12.10f        %12.10f        %12.10f    \n',iter,relchange,abschange,ssq);
				disp(' ')
				disp(' Iterations terminated based on relative change in fit error')
			end
		elseif abschange < tol(2)
			flag = 1;
			if DumpToScreen
				disp(' '),disp('    Iteration    Rel. Change         Abs. Change         sum-sq residuals'),disp(' ')
				fprintf(' %9.0f       %12.10f        %12.10f        %12.10f    \n',iter,relchange,abschange,ssq);
				disp(' ')
				disp(' Iterations terminated based on absolute change in fit error')
			end
		elseif iter > tol(3)-1
			flag = 1;
			if DumpToScreen
				disp(' '),disp('    Iteration    Rel. Change         Abs. Change         sum-sq residuals'),disp(' ')
				fprintf(' %9.0f       %12.10f        %12.10f        %12.10f    \n',iter,relchange,abschange,ssq);
				disp(' ')
				disp(' Iterations terminated based on maximum iterations')
			end
		end
	end
	if rem(iter,Show) == 0&DumpToScreen
		if iter == Show|rem(iter,Show*30) == 0
			disp(' '),disp('    Iteration    Rel. Change         Abs. Change         sum-sq residuals'),disp(' ')
		end
		fprintf(' %9.0f       %12.10f        %12.10f        %12.10f    \n',iter,relchange,abschange,ssq);
	end
	oldssq = ssq;
	
	% Every fifth iteration do a line search if ls == 1
	if (iter/5 == round(iter/5) & ls == 1)
		% Determine the search direction as the difference between the last two estimates
		for ij = 1:order
			searchdir{ij} = x0{ij} - oldx0{ij};
		end
		if iteroffsets
			for ij=1:length(iteroff)
				searchdirOffs{ij} = Offset{ij}.Offsets - OldOffset{ij}.Offsets;
			end
		end      
		% Initialize other variables required for line search
		testmod = x0;
		sflag = 0; 
		i = 0; 
		sd = zeros(10,1); 
		sd(1) = ssq;
		xl = zeros(10,1);
		while sflag == 0
			for k = 1:order
				testmod{k} = testmod{k} + (2^i)*searchdir{k};
			end
			if iteroffsets
				testoffsets = Offset;
				for j=1:length(Offset)
					testoffsets{j}.Offsets = testoffsets{j}.Offsets + (2^i)*searchdirOffs{j};
				end
			end
			% Calculate the fit error on the new test model
			xest = outerm(testmod);
			%Iterative preproc
			if iteroffsets
				XTrilin = xest; % TO SUBTRACT FROM DATA 
				for j=1:length(iteroff) % Add offsets one at a time
					xest = xest + testoffsets{j}.Offsets;
				end
			end
			if DoWeight
				xsq = ((x-outerm(testmod)).*weights).^2;
			else
				xsq = (x-outerm(testmod)).^2;
			end
			if Missing
				xsq(MissId)=0;  
			end
			% Save the difference and the distance along the search direction
			sd(i+2) = sum(xsq(:));
			xl(i+2) = xl(i+1) + 2^i;
			i = i+1;
			% Check to see if a minimum has been exceeded once two new points are calculated
			if i > 1 
				if sd(i+1) > sd(i)
					sflag = 1;
					% Estimate the minimum along the search direction
					xstar = sum((xl([i i+1 i-1]).^2 - xl([i+1 i-1 i]).^2).*sd(i-1:i+1));
					xstar = xstar/(2*sum((xl([i i+1 i-1]) - xl([i+1 i-1 i])).*sd(i-1:i+1)));
					% Save the old model and update the new one
					oldx0 = x0;
					for k = 1:order
						x0{k} = x0{k} + xstar*searchdir{k};
					end
					if iteroffsets
						OldOffset = Offset;
						for j=1:length(iteroff)
							Offset{j}.Offsets = Offset{j}.Offsets + + xstar*searchdirOffs{j};
						end
					end      
					
				end
			end
		end 
		% Calculate the estimate of the input array based on current loads
		xest = outerm(x0);
		%Iterative preproc
		if iteroffsets
			XTrilin = xest; % TO SUBTRACT FROM DATA 
			for j=1:length(iteroff) % Add offsets one at a time
				xest = xest + Offset{j}.Offsets;
			end      
		end
		xsq = xest;
		if DoWeight
			xsq = ((x-xsq).*weights).^2;
		else
			xsq = (x-xsq).^2;
		end
		if Missing
			xsq(MissId)=0;  
		end
		oldssq = sum(xsq(:));
		if Missing 
			x(MissId)=xest(MissId);
			for j = 1:order
				xuf{j} = unfoldmw(x,j)';
			end
		end
		% disp(sprintf('SSQ at xstar is %g',oldssq))
	else
		% Save the last estimates of the loads
		oldx0 = x0;
	end
	% Exchange missing with model estimates
end


% Plot the loadings
if plots ~= 0
	h1 = figure('position',[170 130 512 384],'name','PARAFAC Loadings');
	for i = 1:order
		subplot(order,1,i)
		if ~isempty(scl{i})
			if xsize(i) <= 50
				plot(scl{i},x0{i},'-+')
			else
				plot(scl{i},x0{i},'-')
			end
		else
			if xsize(i) <= 50
				plot(x0{i},'-+')
			else
				plot(x0{i},'-')
			end
		end
		ylabel(sprintf('Dimension %g',i))
		if i == 1
			title('Loadings for Each Dimension')
		end
	end
end

% Calculate and plot the residuals  
dif = (x-outerm(x0)).^2;
res = cell(1,order);
for i = 1:order
	x = dif;
	for j = 1:order
		if i ~= j
			x = sum(x,j);
		end
	end
	x = squeeze(x);
	res{i} = x(:);
	if plots ~= 0
		if i == 1  
			figure('position',[145 166 512 384],'name','PARAFAC Residuals')
		end
		subplot(order,1,i)
		if ~isempty(scl{i})
			if xsize(i) <= 50
				plot(scl{i},res{i},'-+')
			else
				plot(scl{i},res{i},'-')
			end
		else
			if xsize(i) <= 50
				plot(res{i},'-+')
			else
				plot(res{i},'-')
			end
		end
		ylabel(sprintf('Dimension %g',i))
		if i == 1
			title('Residuals for Each Dimension')
		end
	end
end

% Bring the loads back to the front
if plots ~= 0
	figure(h1)
end

% Save the model as a structured array    
model = struct('xname',inputname(1),'name','PARAFAC','date',date,'time',clock,...
	'size',xsize,'nocomp',nocomp,'tol',tol,'final',[relchange abschange iter]);
model.ssq = [tssq ssq];
model.loads = x0;
model.res = res;
model.scl = scl;
%Iterative preproc
if iteroffsets
	for i=1:length(iteroff) % Add offsets one at a time
		model.Offsets.Parameters{i} = Offset{i}.OffsetsReduced;
		model.Offsets.Model{i} = Offset{i}.Offsets;
		model.Offsets.Def{i} = iteroff{i};
	end      
end

function B=unimodal(X,Y,Bold)

% Solves the problem min|Y-XB'| subject to the columns of 
% B are unimodal and nonnegative. The algorithm is iterative
% If an estimate of B (Bold) is given only one iteration is given, hence
% the solution is only improving not least squares
% If Bold is not given the least squares solution is estimated
%
% Copyright 1997
%
% Rasmus Bro
% Royal Veterinary & Agricultural University
% Denmark
% rb@kvl.dk
%
% Reference
% Bro and Sidiropoulos, "Journal of Chemometrics", 1998, 12, 223-247. 


if nargin==3
	B=Bold;
	F=size(B,2);
	for f=1:F
		y=Y-X(:,[1:f-1 f+1:F])*B(:,[1:f-1 f+1:F])';
		beta=pinv(X(:,f))*y;
		B(:,f)=ulsr(beta',1);
	end
else
	F=size(X,2);
	maxit=100;
	B=randn(size(Y,2),F);
	Bold=2*B;
	it=0;
	while norm(Bold-B)/norm(B)>1e-5&it<maxit
		Bold=B;
		it=it+1;
		for f=1:F
			y=Y-X(:,[1:f-1 f+1:F])*B(:,[1:f-1 f+1:F])';
			beta=pinv(X(:,f))*y;
			B(:,f)=ulsr(beta',1);
		end
	end
	if it==maxit
		disp([' UNIMODAL did not converge in ',num2str(maxit),' iterations']);
	end
end


function [b,All,MaxML]=ulsr(x,NonNeg);

% ------INPUT------
%
% x          is the vector to be approximated
% NonNeg     If NonNeg is one, nonnegativity is imposed
%
%
%
% ------OUTPUT-----
%
% b 	     is the best ULSR vector
% All 	     is containing in its i'th column the ULSRFIX solution for mode
% 	     location at the i'th element. The ULSR solution given in All
%            is found disregarding the i'th element and hence NOT optimal
% MaxML      is the optimal (leftmost) mode location (i.e. position of maximum)
%
% ___________________________________________________________
%
%
%               Copyright 1997
%
% Nikos Sidiroupolos
% University of Maryland
% Maryland, US
%
%       &
%
% Rasmus Bro
% Royal Veterinary & Agricultural University
% Denmark
%
% 
% ___________________________________________________________


% This file uses MONREG.M

x=x(:);
I=length(x);
xmin=min(x);
if xmin<0
	x=x-xmin;
end


% THE SUBSEQUENT 
% CALCULATES BEST BY TWO MONOTONIC REGRESSIONS

% B1(1:i,i) contains the monontonic increasing regr. on x(1:i)
[b1,out,B1]=monreg(x);

% BI is the opposite of B1. Hence BI(i:I,i) holds the monotonic
% decreasing regression on x(i:I)
[bI,out,BI]=monreg(flipud(x));
BI=flipud(fliplr(BI));

% Together B1 and BI can be concatenated to give the solution to
% problem ULSR for any modloc position AS long as we do not pay
% attention to the element of x at this position


All=zeros(I,I+2);
All(1:I,3:I+2)=B1;
All(1:I,1:I)=All(1:I,1:I)+BI;
All=All(:,2:I+1);
Allmin=All;
Allmax=All;
% All(:,i) holds the ULSR solution for modloc = i, disregarding x(i),


iii=find(x>=max(All)');
b=All(:,iii(1));
b(iii(1))=x(iii(1));
Bestfit=sum((b-x).^2);
MaxML=iii(1);
for ii=2:length(iii)
	this=All(:,iii(ii));
	this(iii(ii))=x(iii(ii));
	thisfit=sum((this-x).^2);
	if thisfit<Bestfit
		b=this;
		Bestfit=thisfit;
		MaxML=iii(ii);
	end
end

if xmin<0
	b=b+xmin;
end


% Impose nonnegativity
if NonNeg==1
	if any(b<0)
		id=find(b<0);
		% Note that changing the negative values to zero does not affect the
		% solution with respect to nonnegative parameters and position of the
		% maximum.
		b(id)=zeros(size(id))+0;
	end
end



function [b,B,AllBs]=monreg(x);

% Monotonic regression according
% to J. B. Kruskal 64
%
% b     = min|x-b| subject to monotonic increase
% B     = b, but condensed
% AllBs = All monotonic regressions, i.e. AllBs(1:i,i) is the 
%         monotonic regression of x(1:i)
%
%
% Copyright 1997
%
% Rasmus Bro
% Royal Veterinary & Agricultural University
% Denmark
% rb@kvl.dk
%


I=length(x);
if size(x,2)==2
	B=x;
else
	B=[x(:) ones(I,1)];
end

AllBs=zeros(I,I);
AllBs(1,1)=x(1);
i=1;
while i<size(B,1)
	if B(i,1)>B(min(I,i+1),1)
		summ=B(i,2)+B(i+1,2);
		B=[B(1:i-1,:);[(B(i,1)*B(i,2)+B(i+1,1)*B(i+1,2))/(summ) summ];B(i+2:size(B,1),:)];
		OK=1;
		while OK
			if B(i,1)<B(max(1,i-1),1)
				summ=B(i,2)+B(i-1,2);
				B=[B(1:i-2,:);[(B(i,1)*B(i,2)+B(i-1,1)*B(i-1,2))/(summ) summ];B(i+1:size(B,1),:)];
				i=max(1,i-1);
			else
				OK=0;
			end
		end
		bInterim=[];
		for i2=1:i
			bInterim=[bInterim;zeros(B(i2,2),1)+B(i2,1)];
		end
		No=sum(B(1:i,2));
		AllBs(1:No,No)=bInterim;
	else
		i=i+1;
		bInterim=[];
		for i2=1:i
			bInterim=[bInterim;zeros(B(i2,2),1)+B(i2,1)];
		end
		No=sum(B(1:i,2));
		AllBs(1:No,No)=bInterim;
	end
end

b=[];
for i=1:size(B,1)
	b=[b;zeros(B(i,2),1)+B(i,1)];
end


function [x,w] = CrossProdFastNnls(XtX,Xty,tol)
%NNLS	Non-negative least-squares.
%	b = CrossProdFastNnls(XtX,Xty) returns the vector b that solves X*b = y
%	in a least squares sense, subject to b >= 0, given the inputs
%       XtX = X'*X and Xty = X'*y.
%
%	[b,w] = fastnnls(XtX,Xty) also returns dual vector w where
%	w(i) < 0 where b(i) = 0 and w(i) = 0 where b(i) > 0.
%	L. Shure 5-8-87 Copyright (c) 1984-94 by The MathWorks, Inc.
%
%  Revised by:
%	Copyright
%	Rasmus Bro 1995
%	Denmark
%	E-mail rb@kvl.dk
%  According to Bro & de Jong, J. Chemom, 1997, 11, 393-401

% initialize variables


if nargin < 3
	tol = 10*eps*norm(XtX,1)*max(size(XtX));
end
[m,n] = size(XtX);
P = zeros(1,n);
Z = 1:n;
x = P';
ZZ=Z;
w = Xty-XtX*x;

% set up iteration criterion
iter = 0;
itmax = 30*n;

% outer loop to put variables into set to hold positive coefficients
while any(Z) & any(w(ZZ) > tol)
	[wt,t] = max(w(ZZ));
	t = ZZ(t);
	P(1,t) = t;
	Z(t) = 0;
	PP = find(P);
	ZZ = find(Z);
	nzz = size(ZZ);
	z(PP')=(Xty(PP)'/XtX(PP,PP)');
	z(ZZ) = zeros(nzz(2),nzz(1))';
	z=z(:);
	% inner loop to remove elements from the positive set which no longer belong
	
	while any((z(PP) <= tol)) & iter < itmax
		
		iter = iter + 1;
		QQ = find((z <= tol) & P');
		alpha = min(x(QQ)./(x(QQ) - z(QQ)));
		x = x + alpha*(z - x);
		ij = find(abs(x) < tol & P' ~= 0);
		Z(ij)=ij';
		P(ij)=zeros(1,max(size(ij)));
		PP = find(P);
		ZZ = find(Z);
		nzz = size(ZZ);
		z(PP)=(Xty(PP)'/XtX(PP,PP)');
		z(ZZ) = zeros(nzz(2),nzz(1));
		z=z(:);
	end
	x = z;
	w = Xty-XtX*x;
end

x=x(:);





function loads=atld(X,F,Show);



if nargin<3
	Show = 0;
end

xsize = size(X);
order = ndims(X);
RealData = all(isreal(X(:)));

% Initialize with random numbers
for i = 1:order;
	if xsize(i)>=F
		loads{i} = orth(rand(xsize(i),F));
	else
		loads{i} = rand(xsize(i),F);
	end
	invloads{i} = pinv(loads{i});
end

model = outerm(loads);
fit=sum(abs((X(:)-model(:)).^2));
oldfit=2*fit;
maxit=30;
crit=1e-6;
it=0;

while abs(fit-oldfit)/oldfit>crit&it<maxit
	it=it+1;
	oldfit=fit;
   
   
	% Normalize loadings   
	for mode = 2:order
		scale = sqrt(abs(sum(loads{mode}.^2)));
		loads{1}    = loads{1}*diag(scale);
		loads{mode} = loads{mode}*diag(scale.^-1);
		% Scale occasionally to mainly positiviy
		if RealData & (it==1|it==2|rem(it,1)==0)
			scale = sign(sum(loads{mode}));
			loads{1}    = loads{1}*diag(scale);
         loads{mode} = loads{mode}*diag(scale);
         if any(isnan(loads{mode}))
            loads{mode} = rand(size(loads{mode}));
         end
			invloads{mode} = pinv(loads{mode});
		end
	end
	
	if it == 1
		delta = linspace(1,F^(order-1),F);
	end
	
	% Compute new loadings for all modes
	for mode = 1:order;
		
		xprod = X;
		% Multiply pseudo-inverse loadings in all but the mode being estimated
		if mode == 1
			for mulmode = 2:order
				xprod = ntimes(xprod,invloads{mulmode},2,2);
			end
		else
			for j = 1:mode-1
				xprod = ntimes(xprod,invloads{j},1,2);
			end
			for j = mode+1:order
				xprod = ntimes(xprod,invloads{j},2,2);
			end
		end
		
		% Extract first mode loadings from product
		loads{mode} = xprod(:,delta);
		invloads{mode} = pinv(loads{mode});
	end
	
	
	model = outerm(loads);
	fit=sum((X(:)-model(:)).^2);
end
% Normalize loadings   
for mode = 2:order
	scale = sqrt(sum(loads{mode}.^2));
	loads{1}    = loads{1}*diag(scale);
	loads{mode} = loads{mode}*diag(scale.^-1);
end

if Show
	disp(' ')
	disp('    Iteration    sum-sq residuals')
	disp(' ')
	fprintf(' %9.0f       %12.10f    \n',it,fit);
end




function product = ntimes(X,Y,modeX,modeY);


%NTIMES Array multiplication
% X*Y is the array/matrix product of X and Y. These are multiplied across the 
% modeX mode/dimension of X and modeY mode/dimension of Y. The number of levels of X in modeX
% must equal the number of levels in modeY of Y.
%
% The product will be an array of order two less than the sum of the orders of A and B
% and thus works as a straightforward extension of the matrix product. The order
% of the modes in the product are such that the X-modes come first and then the
% Y modes. 
%
% E.g. if X is IxJ and Y is KxL and I equals L, then 
% ntimes(X,Y,1,2) will yield a JxK matrix equivalent to the matrix product
% X'*Y'. If X is IxJxK and Y is LxMxNxO with J = N then 
% ntimes(X,Y,2,3) yields a product of size IxKxLxMxO
% 
%I/O: product = NTIMES(X,Y,modeX,modeY) 
% 
%See also TIMES,MTIMES.

%Copyright Eigenvector Research Inc./Rasmus Bro, 2000
%Rasmus Bro, August 20, 2000

orderX = ndims(X);
orderY = ndims(Y);
xsize  = size(X);
ysize  = size(Y);

X = permute(X,[modeX 1:modeX-1 modeX+1:orderX]);
Y = permute(Y,[modeY 1:modeY-1 modeY+1:orderY]);
xsize2  = size(X);
ysize2  = size(Y);

if size(X,1)~=size(Y,1)
	error(' The number of levels must be the same in the mode across which multiplications is performed')
end

%multiply the matricized arrays
product = reshape(X,xsize2(1),prod(xsize2(2:end)))'*reshape(Y,ysize2(1),prod(ysize2(2:end)));
%reshape to three-way structure
product = reshape(product,[xsize2(2:end) ysize2(2:end)]);



function mwauf = unfoldmw(mwa,order,meth)
%UNFOLDMW Unfolds multiway arrays along specified order
% The inputs are the multiway array to be unfolded (mwa),
% and the dimension number along which to perform the
% unfolding (order). The output is the unfolded array (mwauf).
% This function is used in the development of PARAFAC models
% in the alternating least squares steps.
%
%I/O: mwauf = unfoldmw(mwa,order);
%
%See also: MPCA, OUTER, OUTERM, PARAFAC, TLD

%Copyright Eigenvector Research, Inc. 1998
%bmw
%rb August 2000, speeded up by factor 20


mwasize = size(mwa);
ord = length(mwasize);
mwauf=permute(mwa,[order 1:order-1 order+1:ord]);
mwauf=reshape(mwauf,mwasize(order),prod(mwasize([1:order-1 order+1:ord])));

%Old version
%mwasize = size(mwa);
%ms = mwasize(order);
%po = prod(mwasize);
%ns = po/ms;
%if order ~= 1
%   pod = prod(mwasize(1:order-1));
%end
%mwauf = zeros(ms,ns);
%for i = 1:ms
%   if order == 1
%      mwauf(i,:) = mwa(i:ms:po);
%   else
%      inds = zeros(1,ns); k = 1; fi = (i-1)*pod + 1;
%      for j = 1:ns/pod
%         inds(k:k+pod-1) = fi:fi+pod-1;
%         fi = fi + ms*pod;
%         k = k + pod;
%      end
%      mwauf(i,:) = mwa(inds);
%   end
%end%
%
%end

function mwa = outerm(facts,lo,vect)
%OUTERM Outer product of any number of vectors with multiple factors
%  The input to outerm is a 1 by n cell array (facts), where each cell
%  contains the factors for one of the ways, or orders, with each
%  of the factors being a column in the matrix. Optional inputs
%  are the number of an order to leave out (lo) in the formation
%  of the product, and a flag (vect) which causes the function
%  to not sum and reshape the final factors when set to 1. (This option
%  is used in the alternating least squares steps in PARAFAC.) 
%  The output is the multiway array resulting from multiplying the
%  factors together(mwa), or the strung out individual factors.
%
%I/O: mwa = outerm(facts,lo,vect);
%
%See also: OUTER, PARAFAC, TLD

%Copyright Eigenvector Research, Inc. 1998
%bmw

if nargin < 2
  lo = 0;
end
if nargin < 3
  vect = 0;
end
order = length(facts);
if lo == 0
  mwasize = zeros(1,order);
else
  mwasize = zeros(1,order-1);
end
k = 0;
for i = 1:order
  [m,n] = size(facts{i});
  if i ~= lo
    k = k + 1;
    mwasize(k) = m;
  end
  if i > 1
    if nofac ~= n
	  error('All orders must have the same number of factors')
	end
  else
    nofac = n;
  end
end
mwa = zeros(prod(mwasize),nofac);

for j = 1:nofac
  if lo ~= 1
    mwvect = facts{1}(:,j);
    for i = 2:order
	  if lo ~= i
        %mwvect = kron(facts{i}(:,j),mwvect);
		mwvect = mwvect*facts{i}(:,j)';
		mwvect = mwvect(:);
	  end
    end
  elseif lo == 1
    mwvect = facts{2}(:,j);
	for i = 3:order
      %mwvect = kron(facts{i}(:,j),mwvect);
	  mwvect = mwvect*facts{i}(:,j)';
	  mwvect = mwvect(:);
	end
  end
  mwa(:,j) = mwvect;
end
% If vect isn't one, sum up the results of the factors and reshape
if vect ~= 1
  mwa = sum(mwa,2);
  mwa = reshape(mwa,mwasize);
end

function model = tld(x,ncomp,scl,plots)
%TLD Trilinear decomposition.
%  The Trilinear decomposition can be used to decompose
%  a 3-way array as the summation over the outer product
%  of triads of vectors. The inputs are the 3 way array
%  (x) and the number of components to estimate (ncomp),
%  Optional input variables include a 1 by 3 cell array 
%  containing scales for plotting the profiles in each
%  order (scl) and a flag which supresses the plots when
%  set to zero (plots). The output of TLD is a structured
%  array (model) containing all of the model elements
%  as follows:
%
%     xname: name of the original workspace input variable
%      name: type of model, always 'TLD'
%      date: model creation date stamp
%      time: model creation time stamp
%      size: size of the original input array
%    nocomp: number of components estimated
%     loads: 1 by 3 cell array of the loadings in each dimension
%       res: 1 by 3 cell array residuals summed over each dimension
%       scl: 1 by 3 cell array with scales for plotting loads
%
%  Note that the model loadings are presented as unit vectors
%  for the first two dimensions, remaining scale information is
%  incorporated into the final (third) dimension. 
%
%I/O: model = tld(x,ncomp,scl,plots);
%
%See also: GRAM, MWFIT, OUTER, OUTERM, PARAFAC

%Copyright Eigenvector Research, Inc. 1998
%By Barry M. Wise
%Modified April, 1998 BMW
%Modified May, 2000 BMW


dx = size(x);

[min_dim,min_mode] = min(dx);
shift_mode = min_mode;
if shift_mode == 3
   shift_mode = 0;
end
x = shiftdim(x,shift_mode);
dx = size(x);


if (nargin < 3 | ~strcmp(class(scl),'cell'))
  scl = cell(1,3);
end
if nargin < 4
  plots = 1;
end
xu = reshape(x,dx(1),dx(2)*dx(3));
if dx(1) > dx(2)*dx(3)
  [u,s,v] = svd(xu,0);
else
  [v,s,u] = svd(xu',0);
end
uu = u(:,1:ncomp);
xu = zeros(dx(2),dx(1)*dx(3));
for i = 1:dx(1)
  xu(:,(i-1)*dx(3)+1:i*dx(3)) = squeeze(x(i,:,:));
end
if dx(2) > dx(1)*dx(3)
  [u,s,v] = svd(xu,0);
else
  [v,s,u] = svd(xu',0);
end
vv = u(:,1:ncomp);

xu = zeros(dx(3),dx(1)*dx(2));
for i = 1:dx(2)
  xu(:,(i-1)*dx(1)+1:i*dx(1)) = squeeze(x(:,i,:))';
end
if dx(3) > dx(1)*dx(2)
  [u,s,v] = svd(xu,0);
else
  [v,s,u] = svd(xu',0);
end
ww = u(:,1:2);
clear u s v

g1 = zeros(ncomp,ncomp,dx(3));
for i = 1:dx(3)
  g1(:,:,i) = uu'*squeeze(x(:,:,i))*vv;
end
g2 = g1;
for i = 1:dx(3);
  g1(:,:,i) = g1(:,:,i)*ww(i,2);
  g2(:,:,i) = g2(:,:,i)*ww(i,1);
end
g1 = sum(g1,3);
g2 = sum(g2,3);
[aa,bb,qq,zz,ev] = qz(g1,g2);
if ~isreal(ev)
  disp('  ')
  disp('Imaginary solution detected')
  disp('Rotating Eigenvectors to nearest real solution')
  ev=simtrans(aa,bb,ev);
end
ord1 = uu*(g1)*ev;
ord2 = vv*pinv(ev');
norms1 = sqrt(sum(ord1.^2));
norms2 = sqrt(sum(ord2.^2));
ord1 = ord1*inv(diag(norms1));
ord2 = ord2*inv(diag(norms2));
sf1 = sign(mean(ord1));
if any(sf1==0)
  sf1(find(sf1==0)) = 1;
end
ord1 = ord1*diag(sf1);
sf2 = sign(mean(ord2));
if any(sf2==0)
  sf2(find(sf2==0)) = 1;
end
ord2 = ord2*diag(sf2);
ord3 = zeros(dx(3),ncomp);
xu = zeros(dx(1)*dx(2),ncomp);
for i = 1:ncomp
  xy = ord1(:,i)*ord2(:,i)';
  xu(:,i) = xy(:);
end
for i = 1:dx(3)
  y = squeeze(x(:,:,i));
  ord3(i,:) = (xu\y(:))';
end

if shift_mode
   if shift_mode==1
      ord4 = ord1;
      ord1 = ord3;
      ord3 = ord2;
      ord2 = ord4;
   else
      ord4 = ord1;
      ord1 = ord2;
      ord2 = ord3;
      ord3 = ord4;
   end
   x = shiftdim(x,3-shift_mode);
   dx = size(x);
end

if plots ~= 0
  h1 = figure('position',[170 130 512 384],'name','TLD Loadings');
  subplot(3,1,1)
  if ~isempty(scl{1})
    if dx(1) < 50
      plot(scl{1},ord1,'-+')
    else
      plot(scl{1},ord1,'-')
    end
  else
    if dx(1) < 50
      plot(ord1,'-+')
    else
      plot(ord1,'-')
    end
  end
  title('Profiles in First Order')
  subplot(3,1,2)
  if ~isempty(scl{2})
    if dx(2) < 50
      plot(scl{2},ord2,'-+')
    else
      plot(scl{2},ord2,'-')
    end
  else
    if dx(2) < 50
      plot(ord2,'-+')
    else
      plot(ord2,'-')
    end
  end
  title('Profiles in Second Order')
  subplot(3,1,3)
  if ~isempty(scl{3})
    if dx(3) < 50
      plot(scl{3},ord3,'-+')
    else
      plot(scl{3},ord3,'-')
    end
  else
    if dx(3) < 50
      plot(ord3,'-+')
    else
      plot(ord3,'-')
    end
  end
  title('Profiles in Third Order')
end
loads = cell(1,3);
loads{1} = ord1;
loads{2} = ord2;
loads{3} = ord3;
xhat = outerm(loads);
dif = (x-xhat).^2;
res = cell(1,3);
res{1} = sum(sum(dif,3),2);
res{2} = sum(sum(dif,3),1)';
res{3} = squeeze(sum(sum(dif,2),1));
if plots ~= 0
  figure('position',[145 166 512 384],'name','TLD Residuals')
  subplot(3,1,1)
  if ~isempty(scl{1})
    if dx(1) < 50
      plot(scl{1},res{1},'-+')
    else
      plot(scl{1},res{1},'-')
    end
  else
    if dx(1) < 50
      plot(res{1},'-+')
    else
      plot(res{1},'-')
    end
  end
  title('Residuals in First Order')
  subplot(3,1,2)
  if ~isempty(scl{2})
    if dx(2) < 50
      plot(scl{2},res{2},'-+')
    else
      plot(scl{2},res{2},'-')
    end
  else
    if dx(2) < 50
      plot(res{2},'-+')
    else
      plot(res{2},'-')
    end
  end
  title('Residuals in Second Order')
  subplot(3,1,3)
  if ~isempty(scl{3})
    if dx(3) < 50
      plot(scl{3},res{3},'-+')
    else
      plot(scl{3},res{3},'-')
    end
  else
    if dx(3) < 50
      plot(res{3},'-+')
    else
      plot(res{3},'-')
    end
  end
  title('Residuals in Third Order')
end

% Bring the loads back to the front
if plots ~= 0
  figure(h1)
end

model = struct('xname',inputname(1),'name','TLD','date',date,'time',clock,...
  'size',dx,'nocomp',ncomp);
model.loads = loads;
model.res = res;
model.scale = scl;


function Vdd=simtrans(aa,bb,ev);
%SIMTRANS Similarity transform to rotate eigenvectors to real solution
Lambda = diag(aa)./diag(bb);
n=length(Lambda);
[t,o]=sort(Lambda);
Lambda(n:-1:1)=Lambda(o);
ev(:,n:-1:1)=ev(:,o);

Theta = angle(ev);
Tdd = zeros(n);
Td = zeros(n);
ii = sqrt(-1);

k=1;
while k <= n
  if k == n
    Tdd(k,k)=1;
    Td(k,k)=(exp(ii*Theta(k,k)));
    k = k+1;
  elseif abs(Lambda(k))-abs(Lambda(k+1)) > (1e-10)*abs(Lambda(k)) 
    %Not a Conjugate Pair
    Tdd(k,k)=1;
    Td(k,k)=(exp(ii*Theta(k,k)));
    k = k+1;
  else 
    %Is a Conjugate Pair
    Tdd(k:k+1,k:k+1)=[1, 1; ii, -ii];
    Td(k,k)=(exp(ii*0));  
    Td(k+1,k+1)=(exp(ii*(Theta(k,k+1)+Theta(k,k))));
    k = k+2;
  end
end
Vd = ev*pinv(Td);
Vdd = Vd*pinv(Tdd);
if imag(Vdd) < 1e-3
   Vdd = real(Vdd);
end
