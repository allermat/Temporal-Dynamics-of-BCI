function fakeData = generateFakeResponses(cfg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Parsing input
p = inputParser;

validModels = {'bci','fus','segA','segV','taskRel'};
validResp = {'cont','discr'};
validDecisionFun = {'mdlAvg','mdlSel','probMatch'};

addParameter(p,'model','',@(x) any(validatestring(x,validModels)));
addParameter(p,'param',struct,@(x) validateattributes(x,{'struct'},{'nonempty'}));
addParameter(p,'nCond',1,@(x) validateattributes(x,{'numeric'},...
    {'integer','positive'}));
addParameter(p,'nTrialPerCond',1,@(x) validateattributes(x,{'numeric'},...
    {'integer','positive'}));
addParameter(p,'resp','cont',@(x) any(validatestring(x,validResp)));
addParameter(p,'nBin',[],@(x) validateattributes(x,{'numeric'},...
    {'scalar','positive'}));
addParameter(p,'decisionFun','',@(x) any(validatestring(x,validDecisionFun)));

parse(p,cfg);

model = p.Results.model;
param = p.Results.param;
nCond = p.Results.nCond;
nTrialPerCond = p.Results.nTrialPerCond;
resp = p.Results.resp;
nBin = p.Results.nBin;
decisionFun = p.Results.decisionFun;

rng('shuffle');
s = rng;

switch model
    
    case 'bci'
        
        % Original paramters
        p_common = param.p_common;
        sigP = param.sigP;
        sigA = param.sigA;
        sigV = param.sigV;
        if strcmp(resp,'cont'), kW = param.kW; end
        if isempty(decisionFun); decisionFun = 'mdlAvg'; end
        muP = param.muP;
        
        isCommonSource = binornd(1,p_common,nCond,1);
        [sA,sV] = deal(NaN(nCond,1));
        for i = 1:numel(isCommonSource)
            if isCommonSource(i)
                % Drawing one common stimulus location (condition) from 
                % normal distribution
                sA(i) = (randn*sigP)+muP;
                sV(i) = sA(i);
            else
                % Drawing two independent stimulus locations (condition) from 
                % normal distribution
                sA(i) = (randn*sigP)+muP;
                sV(i) = (randn*sigP)+muP;
            end
        end
        
        % Replicating each condition nTiralPerCond times
        sA = repmat(sA,nTrialPerCond,1);
        sV = repmat(sV,nTrialPerCond,1);
        
        [sV_hat,sA_hat,xV,xA] = deal(NaN(size(sA)));
        
        % Variances of A and V and prior
        varV = sigV^2;
        varA = sigA^2;
        varP = sigP^2;
        
        % Variances of estimates given common or independent causes
        varVA_hat = 1/(1/varV + 1/varA + 1/varP);
        varV_hat = 1/(1/varV + 1/varP);
        varA_hat = 1/(1/varA + 1/varP);
        
        % Variances used in computing probability of common or independent causes
        var_common = varV * varA + varV * varP + varA * varP;
        varV_indep = varV + varP;
        varA_indep = varA + varP;
        
        for i = 1:size(sA,1)
            % Generating fake data by adding random noise to the stimulus
            % locations
            xV(i) = sV(i) + sigV * randn;
            xA(i) = sA(i) + sigA * randn;
            
            % Estimates given common source
            s_hat_common = (xV(i)/varV + xA(i)/varA + muP/varP) * varVA_hat;
            sV_hat_indep = (xV(i)/varV + muP/varP) * varV_hat;
            sA_hat_indep = (xA(i)/varA + muP/varP) * varA_hat;
            
            % Probability of common or independent causes
            quad_common = (xV(i)-xA(i)).^2 * varP + (xV(i)-muP).^2 * varA + (xA(i)-muP).^2 * varV;
            quadV_indep = (xV(i)-muP).^2;
            quadA_indep = (xA(i)-muP).^2;
            
            % Likelihood of observations (xV,xA) given C, for C=1 and C=2
            likelihood_common = exp(-quad_common/(2*var_common))/(2*pi*sqrt(var_common));
            likelihoodV_indep = exp(-quadV_indep/(2*varV_indep))/sqrt(2*pi*varV_indep);
            likelihoodA_indep = exp(-quadA_indep/(2*varA_indep))/sqrt(2*pi*varA_indep);
            likelihood_indep =  likelihoodV_indep .* likelihoodA_indep;
            
            % Posterior probability of C given observations (xV,xA)
            post_common = likelihood_common * p_common;
            post_indep = likelihood_indep * (1-p_common);
            pC = post_common./(post_common + post_indep);
            
            switch decisionFun
                
                case 'mdlAvg'
                    % Mean of posterior - Model averaging
                    % Overall estimates: weighted averages
                    sV_hat(i) = pC .* s_hat_common + (1-pC) .* sV_hat_indep;
                    sA_hat(i) = pC .* s_hat_common + (1-pC) .* sA_hat_indep;
                    
                case 'mdlSel'
                    % Model selection instead of model averaging
                    % averaging
                    sV_hat(i) = (pC>0.5).*s_hat_common + (pC<=0.5).*sV_hat_indep;
                    sA_hat(i) = (pC>0.5).*s_hat_common + (pC<=0.5).*sA_hat_indep;
                    
                case 'probMatch'
                    % Probability matching
                    thresh = rand(1,1); % to enable probability matching
                    sV_hat(i) = (pC>thresh).*s_hat_common + (pC<=thresh).*sV_hat_indep;
                    sA_hat(i) = (pC>thresh).*s_hat_common + (pC<=thresh).*sA_hat_indep;
                    
            end
            
        end
        
        if strcmp(resp,'cont')
            % Adding response noise according to the kernel width
            respA = sA_hat + (randn(size(sA_hat))*kW);
            respV = sV_hat + (randn(size(sV_hat))*kW);
        else
            % Discretizing the fake responses
            [~,respLoc] = hist(cat(1,sA_hat,sV_hat),nBin);
            [~,tV] = min(abs(repmat(sV_hat,1,length(respLoc)) - repmat(respLoc,nCond*nTrialPerCond,1)),[],2);
            [~,tA] = min(abs(repmat(sA_hat,1,length(respLoc)) - repmat(respLoc,nCond*nTrialPerCond,1)),[],2);
            respV = respLoc(tV)';
            respA = respLoc(tA)';
        end
        
        fakeData = table(sV,sA,respV,respA,'VariableNames',{'locV','locA','respV','respA'});
        
    case 'fus'
        
        % Original paramters
        sigP = param.sigP;
        sigA = param.sigA;
        sigV = param.sigV;
        if strcmp(resp,'cont'), kW = param.kW; end
        muP = param.muP;
        
        % Drawing stimulus locations (conditions) from normal distribution
        sA = (randn(nCond,1)*sigP)+muP;
        % Replicating each condition nTiralPerCond times
        sA = repmat(sA,1,nTrialPerCond);
        sA = sA(:);
        sV = sA;
        [s_hat_common,xV,xA] = deal(NaN(size(sA)));
        
        % Variances of A and V and prior
        varV = sigV^2;
        varA = sigA^2;
        varP = sigP^2;
        
        % Variances of estimates given common or independent causes
        varVA_hat = 1/(1/varV + 1/varA + 1/varP);
        
        for i = 1:size(sA,1)
            % Generating fake data by adding random noise to the stimulus
            % locations
            xV(i) = sV(i) + sigV * randn;
            xA(i) = sA(i) + sigA * randn;
            
            % Estimates given common source
            s_hat_common(i) = (xV(i)/varV + xA(i)/varA + muP/varP) * varVA_hat;
        end
        
        if strcmp(resp,'cont')
            % Adding response noise according to the kernel width
            respAV = s_hat_common + (randn(size(s_hat_common))*kW);
        else
            % Discretizing the fake responses
            [~,respLoc] = hist(s_hat_common,nBin);
            [~,tV] = min(abs(repmat(s_hat_common,1,length(respLoc)) - repmat(respLoc,nCond*nTrialPerCond,1)),[],2);
            respAV = respLoc(tV);
            respAV = respAV';
        end
        
        fakeData = table(sV,sA,respAV,respAV,'VariableNames',{'locV','locA','respV','respA'});
        
    otherwise
        error('bci:generateFakeResponses:modelNotImplemented',...
            'The requested model is not yet implemented.')
end



end

