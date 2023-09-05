function [B,b,t] = matching_gen_model_mult_old(A,PD,m,modelvar,PDexpo,gam,epsilon)
% GENERATIVE_MODEL          Run generative model code
%
%   B = GENERATIVE_MODEL(A,D,m,modeltype,modelvar,params)
%
%   Generates synthetic networks using the models described in the study by
%   Betzel et al (2016) in Neuroimage.
%
%   Inputs:
%           A,          binary network of seed connections
%           D,          Euclidean distance/fiber length matrix
%           O,          Other distance metrics, with each as the third
%                       dimension in a matrix
%           m,          number of connections that should be present in
%                       final synthetic network
%           modeltype,  specifies the generative rule (see below)
%           modelvar,   specifies whether the generative rules are based on
%                       power-law or exponential relationship
%                       ({'powerlaw'}|{'exponential})
%           eta,        the parameter controlling distance
%           gam,        the parameter controlling topology
%           epsilon,    the baseline probability of forming a particular
%                       connection (should be a very small number
%                       {default = 1e-5}).
%
%   Output:
%           B,          an adjacency matrix
%
%
%   Full list of model types:
%   (each model type realizes a different generative rule)
%
%       1.  'sptl'          spatial model
%       2.  'neighbors'     number of common neighbors
%       3.  'matching'      matching index
%       4.  'clu-avg'       average clustering coeff.
%       5.  'clu-min'       minimum clustering coeff.
%       6.  'clu-max'       maximum clustering coeff.
%       7.  'clu-diff'      difference in clustering coeff.
%       8.  'clu-prod'      product of clustering coeff.
%       9.  'deg-avg'       average degree
%       10. 'deg-min'       minimum degree
%       11. 'deg-max'       maximum degree
%       12. 'deg-diff'      difference in degree
%       13. 'deg-prod'      product of degree
%       14. 'com'           communicability
%
%   Example usage:
%
%       load demo_generative_models_data
%
%       % get number of bi-directional connections
%       m = nnz(A)/2;
% 
%       % get cardinality of network
%       n = length(A);
% 
%       % set model type
%       modeltype = 'neighbors';
% 
%       % set whether the model is based on powerlaw or exponentials
%       modelvar = [{'powerlaw'},{'powerlaw'}];
% 
%       % choose some model parameters
%       params = [-2,0.2; -5,1.2; -1,1.5];
%       nparams = size(params,1);
% 
%       % generate synthetic networks
%       B = generative_model(Aseed,D,m,modeltype,modelvar,params);
%
%       % store them in adjacency matrix format
%       Asynth = zeros(n,n,nparams);
%       for i = 1:nparams; 
%           a = zeros(n); a(B(:,i)) = 1; a = a + a'; 
%           Asynth(:,:,i) = a; 
%       end
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015

if ~exist('epsilon','var')
    epsilon = 1e-6;
end

n = length(A);

% Perform the multiplication of Fo and (D.^eta) or exp(eta*D) as these
% values will not change

nPD = length(PD);

Df = zeros(n,n,nPD);

mv1 = modelvar{1};

if iscell(mv1)
    for ii = 1:nPD
         switch mv1{ii}
            case 'powerlaw'
                    Df(:,:,ii) = PD{ii}.^PDexpo(ii);
            case 'exponential'       
                    Df(:,:,ii) = exp(PDexpo(ii)*(PD{ii}));
         end    
    end    
else
    switch mv1
        case 'powerlaw'
            for i = 1:nPD
                Df(:,:,i) = PD{i}.^PDexpo(i);
            end
        case 'exponential'       
            for i = 1:nPD
                Df(:,:,i) = exp(PDexpo(i)*(PD{i}));  
            end
     end
end

Fd = prod(Df,3);

Kseed = matching_ind_und(A);
[b,t] = fcn_matching(A,Kseed,Fd,m,gam,modelvar,epsilon);
        

B = zeros(n);
B(b) = 1;
B = B + B';

function [b,t] = fcn_matching(A,K,Fd,m,gam,modelvar,epsilon)
K = K + epsilon;
n = length(Fd);
mseed = nnz(A)/2;

mv2 = modelvar{2};

switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);

b = zeros(m,1);
b(1:mseed) = find(A(indx));

t = zeros(length((mseed + 1):m),1);

for ii = (mseed + 1):m
    tic
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(ii) = r;
    uu = u(r);
    vv = v(r);
    
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    
    updateuu = find(A*A(:,uu));
    updateuu(updateuu == uu) = [];
    updateuu(updateuu == vv) = [];
    
    updatevv = find(A*A(:,vv));
    updatevv(updatevv == uu) = [];
    updatevv(updatevv == vv) = [];
    
    c1 = [A(:,uu)', A(uu,:)];
    for i = 1:length(updateuu)
        j = updateuu(i);
        c2 = [A(:,j)' A(j,:)];
        use = ~(~c1&~c2);
        use(uu) = 0;  use(uu+n) = 0;
        use(j) = 0;  use(j+n) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            K(uu,j) = epsilon;
            K(j,uu) = epsilon;
        else
            K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
            K(j,uu) = K(uu,j);
        end
        
    end
    
    c1 = [A(:,vv)', A(vv,:)];
    for i = 1:length(updatevv)
        j = updatevv(i);
        c2 = [A(:,j)' A(j,:)];
        use = ~(~c1&~c2);
        use(vv) = 0;  use(vv+n) = 0;
        use(j) = 0;  use(j+n) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            K(vv,j) = epsilon;
            K(j,vv) = epsilon;
        else
            K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
            K(j,vv) = K(vv,j);
        end
    end
    switch mv2
        case 'powerlaw'
            Fk = K.^gam;
        case 'exponential'
            Fk = exp(gam*K);
    end
    Ff = Fd.*Fk.*~A;
    P = Ff(indx);
    t(ii) = toc;
end
b = indx(b);