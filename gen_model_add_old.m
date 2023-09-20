function [B,b] = gen_model_add_old(A,PDMs,m,modeltype,modelvar,PDMsParam,TopoParam,normType,epsilon)
% gen_model_add_normmax          Run generative model code for the additive
% model normalising by the max for each term
%
%   Generates synthetic networks using the models described in the study by
%   Oldham et al (2022; 10.1126/sciadv.abm6127)
%
%   Inputs:
%           A,          binary network of seed connections
%           PDMs,       Euclidean distance/fiber length/node similarity
%                       matrix. Multiple can be input either as a cell, 
%                       where each cell contains a different matrix or as a
%                       3D matrix (n*n*nPD, where n is the number of nodes
%                       and nPD is the number of PD matrices).
%           m,          number of connections that should be present in
%                       final synthetic network
%           modeltype,  specifies the generative rule (see below)
%           modelvar,   specifies whether the generative rules are based on
%                       power-law or exponential relationship
%                       ({'powerlaw'}|{'exponential})
%           PDMsParam,  the parameters controlling the values in PDMs. 
%                       PDMsParam(1,:) are the exponent values (i.e., the
%                       parameter in the power-law or exponential). The
%                       index of each column corresponds to the index in
%                       PDMs. PDMsParam(2,:) are the alpha values                       
%           TopoParam,  the parameters controlling the topological term. 
%                       TopoParam(1,1) are the exponent values (i.e., the
%                       parameter in the power-law or exponential). The
%                       TopoParam(2,1) are the alpha values 
%           normType,   either 'max' or 'sum', this will be the type of
%                       normalisation performed within each respective
%                       term. The default is 'max'
%           epsilon,    the baseline probability of forming a particular
%                       connection (should be a very small number
%                       {default = 0}).
%
%   Output:
%           B,          an adjacency matrix
%           b,          a vector giving the index of each edge in B. Note
%                       that the ordering of b shows which edges formed
%                       first (e.g., b(1) was the fiorst edge to form, b(2)
%                       the second etc etc).
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
%       How to convert b to B:
%       n = length(A); B = zeros(n); B(b(:,i)) = 1; B = B + B'; 
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%              Oldham et al (2022) Science Advances 10.1126/sciadv.abm6127
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015
%   Edited by Stuart Oldham, Monash University 2021

if ~exist('epsilon','var')
    epsilon = 0;
end

if ~exist('normType','var')
    normType = 'max';
end

n = length(A);

if iscell(PDMs)
    PDMs_ = PDMs;
    nPDMs = length(PDMs);
    PDMs = zeros(n,n,nPDMs);
    for i = 1:nPDMs
        PDMs(:,:,i) = PDMs_{i};        
    end
else
    nPDMs = size(PDMs,3);
end

PDMs=bsxfun(@times, ~eye(size(PDMs(:,:,1))), PDMs);

if nPDMs == 1
    modelvar{1} = {modelvar{1}};
end

if size(PDMsParam,1) == 1
    PDMsParam(2,nPDMs) = 1;
end

if size(TopoParam,1) == 1
    TopoParam(2,1) = 1;
end

PDMsParam(isnan(PDMsParam)) = 0;
TopoParam(isnan(TopoParam)) = 0;

switch modeltype

    case 'clu-avg'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@plus,clu(:,ones(1,n)),clu')/2;
        b = fcn_clu(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'clu-diff'
        clu = clustering_coef_bu(A);
        Kseed = abs(bsxfun(@minus,clu(:,ones(1,n)),clu'));
        b = fcn_clu(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'clu-max'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@max,clu(:,ones(1,n)),clu');
        b = fcn_clu(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'clu-min'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@min,clu(:,ones(1,n)),clu');
        b = fcn_clu(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'clu-prod'
        clu = clustering_coef_bu(A);
        Kseed = clu*clu';
        b = fcn_clu(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'deg-avg'
        kseed = sum(A,2);
        Kseed = bsxfun(@plus,kseed(:,ones(1,n)),kseed')/2;
        b = fcn_deg(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'deg-diff'
        kseed = sum(A,2);
        Kseed = abs(bsxfun(@minus,kseed(:,ones(1,n)),kseed'));
        b = fcn_deg(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'deg-max'
        kseed = sum(A,2);
        Kseed = bsxfun(@max,kseed(:,ones(1,n)),kseed');
        b = fcn_deg(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'deg-min'
        kseed = sum(A,2);
        Kseed = bsxfun(@min,kseed(:,ones(1,n)),kseed');
        b = fcn_deg(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'deg-prod'
        kseed = sum(A,2);
        Kseed = (kseed*kseed').*~eye(n);
        b = fcn_deg(modeltype,A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'neighbors'
        Kseed = (A*A).*~eye(n);
        b = fcn_nghbrs(A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'matching'
        Kseed = matching_ind_und(A);
        b = fcn_matching(A,Kseed,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType);

    case 'sptl'
        b = fcn_sptl(A,PDMs,m,PDMsParam,modelvar,normType);

end

B = zeros(n);
B(b) = 1;
B = B + B';

function b = fcn_clu(modeltype,A,K,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType)
K = K + epsilon;

n = length(A);
mseed = nnz(A)/2;
A = A > 0;

mv2 = modelvar{2};

[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;

PDMsmv = modelvar{1};
PDf = zeros(size(PDMs));
Df1 = zeros(size(PDMs));
for j = 1:length(PDMsmv)
    mv1 = PDMsmv{j};
       
    switch mv1
        case 'powerlaw'
            Df1(:,:,j) = PDMs(:,:,j).^PDMsParam(1,j);
        case 'exponential'       
            Df1(:,:,j) = exp(PDMsParam(1,j)*(PDMs(:,:,j))).*~eye(n);                 
    end
    Df = Df1(:,:,j);
    switch normType
    case 'max'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./max(Df.*~A,[],'all'));
    case 'sum'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./sum(Df.*~A,'all'));        
    end
end

Fd = sum(PDf,3);

switch mv2
    case 'powerlaw'
        Kf = K.^TopoParam(1,1);
    case 'exponential'       
        Kf = exp(TopoParam(1,1)*K);    
end
Kf(isinf(Kf)) = 0;
switch normType
    case 'max'
        Fk = TopoParam(2,1)*(Kf./max(Kf.*~A,[],'all'));
    case 'sum'
        Fk = TopoParam(2,1)*(Kf./sum(Kf.*~A,'all'));
end

Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;

c = clustering_coef_bu(A);
k = sum(A,2);

Ff = (Fd+Fk).*~A;

P = Ff(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    
    b(i) = r;
    EdgesTriu(b(i)) = 1;
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    k([uu,vv]) = k([uu,vv]) + 1;
    bu = A(uu,:);
    su = A(bu,bu);
    bv = A(vv,:);
    sv = A(bv,bv);
    bth = bu & bv;
    c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
    c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
    c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
    c(k <= 1) = 0;
    bth([uu,vv]) = true;
    
    switch modeltype
        case 'clu-avg'
    
            K(:,bth) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')/2 + epsilon;
            K(bth,:) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')'/2 + epsilon;

        case 'clu-diff'

            K(:,bth) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)')) + epsilon;
            K(bth,:) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)'))' + epsilon;

        case 'clu-max'
            
            K(:,bth) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
            K(bth,:) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;
            
        case 'clu-min'
            
            K(:,bth) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
            K(bth,:) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;
    
        case 'clu-prod'
            
            K(bth,:) = (c(bth,:)*c') + epsilon;
            K(:,bth) = (c*c(bth,:)') + epsilon;
    end
    
    switch mv2
        case 'powerlaw'
            Fk = K.^TopoParam(1,1);
        case 'exponential'       
            Fk = exp(TopoParam(1,1)*K);    
    end
    Fk(isinf(Fk)) = 0;
    
    switch normType
        case 'max'
            for j = 1:length(PDMsmv)
                Df = Df1(:,:,j);
                PDf(:,:,j) = PDMsParam(2,j)*(Df./max(Df.*~A,[],'all'));
            end
            Fd = sum(PDf,3);
            MaxNorm = max(Fk.*~A,[],'all');
            if MaxNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = TopoParam(2,1)*(Fk./MaxNorm).*~A;   
            end
        case 'sum'
            for j = 1:length(PDMsmv)
                Df = Df1(:,:,j);
                PDf(:,:,j) = PDMsParam(2,j)*(Df./sum(Df.*~A,'all'));
            end
            SumNorm = sum(Fk.*~A,'all')/2;
            if SumNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = TopoParam(2,1)*(Fk./SumNorm).*~A;   
            end
    end
        
    Ff = (Fd+TopoTerm).*~A;
    P = Ff(indx);
    
end
b = indx(b);


function b = fcn_deg(modeltype,A,K,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType)
n = length(A);
mseed = nnz(A)/2;
k = sum(A,2);
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;

K = K + epsilon;
mv2 = modelvar{2};

PDMsmv = modelvar{1};
PDf = zeros(size(PDMs));
Df1 = zeros(size(PDMs));
for j = 1:length(PDMsmv)
    mv1 = PDMsmv{j};
       
    switch mv1
        case 'powerlaw'
            Df1(:,:,j) = PDMs(:,:,j).^PDMsParam(1,j);
        case 'exponential'       
            Df1(:,:,j) = exp(PDMsParam(1,j)*(PDMs(:,:,j))).*~eye(n);                 
    end
    Df = Df1(:,:,j);
    switch normType
    case 'max'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./max(Df.*~A,[],'all'));
    case 'sum'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./sum(Df.*~A,'all'));        
    end
end

Fd = sum(PDf,3);

switch mv2
    case 'powerlaw'
        Kf = K.^TopoParam(1,1);
    case 'exponential'       
        Kf = exp(TopoParam(1,1)*K);    
end
Kf(isinf(Kf)) = 0;
switch normType
    case 'max'
        Fk = TopoParam(2,1)*(Kf./max(Kf.*~A,[],'all'));
    case 'sum'
        Fk = TopoParam(2,1)*(Kf./sum(Kf.*~A,'all'));
end

Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    
P = (Fd(indx)+Fk(indx)).*~A(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    
    
    switch modeltype
        case 'deg-avg'
    
            K(:,w) = [((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon];
            K(w,:) = ([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon])';
    
        case 'deg-diff'
    
            K(:,w) = (abs([k - k(w(1)), k - k(w(2))]) + epsilon);
            K(w,:) = (abs([k - k(w(1)), k - k(w(2))]) + epsilon)';
    
        case 'deg-min'
            
            K(:,w) = [min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon];
            K(w,:) = [min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon]'; 
            
            
        case 'deg-max'
            
            K(:,w) = [max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon];
            K(w,:) = [max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon]'; 
            
        case 'deg-prod'
            
            K(:,w) = [k*k(w(1)) + epsilon, k*k(w(2)) + epsilon];
            K(w,:) = [k*k(w(1)) + epsilon, k*k(w(2)) + epsilon]';
    
    end
    
    b(i) = r;
    EdgesTriu(b(i)) = 1;
        
    switch mv2
        case 'powerlaw'
            Fk = K.^TopoParam(1,1);
        case 'exponential'       
            Fk = exp(TopoParam(1,1)*K);    
    end
    Fk(isinf(Fk)) = 0;
    
    switch normType
        case 'max'
            for j = 1:length(PDMsmv)
                Df = Df1(:,:,j);
                PDf(:,:,j) = PDMsParam(2,j)*(Df./max(Df.*~A,[],'all'));
            end
            Fd = sum(PDf,3);
            MaxNorm = max(Fk.*~A,[],'all');
            if MaxNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = TopoParam(2,1)*(Fk./MaxNorm).*~A;   
            end
        case 'sum'
            for j = 1:length(PDMsmv)
                Df = Df1(:,:,j);
                PDf(:,:,j) = PDMsParam(2,j)*(Df./sum(Df.*~A,'all'));
            end
            SumNorm = sum(Fk.*~A,'all')/2;
            if SumNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = TopoParam(2,1)*(Fk./SumNorm).*~A;   
            end
    end
        
    Fp = Fd + TopoTerm;

    P = Fp(indx);
    
    P(b(1:i)) = 0;
    
end
b = indx(b);

function b = fcn_nghbrs(A,K,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType)
K = K + epsilon;

n = length(A);
mseed = nnz(A)/2;
A = A > 0;

mv2 = modelvar{2};

[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;

PDMsmv = modelvar{1};
PDf = zeros(size(PDMs));
Df1 = zeros(size(PDMs));
for j = 1:length(PDMsmv)
    mv1 = PDMsmv{j};
       
    switch mv1
        case 'powerlaw'
            Df1(:,:,j) = PDMs(:,:,j).^PDMsParam(1,j);
        case 'exponential'       
            Df1(:,:,j) = exp(PDMsParam(1,j)*(PDMs(:,:,j))).*~eye(n);                 
    end
    Df = Df1(:,:,j);
    switch normType
    case 'max'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./max(Df.*~A,[],'all'));
    case 'sum'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./sum(Df.*~A,'all'));        
    end
end

Fd = sum(PDf,3);

switch mv2
    case 'powerlaw'
        Kf = K.^TopoParam(1,1);
    case 'exponential'       
        Kf = exp(TopoParam(1,1)*K);    
end
Kf(isinf(Kf)) = 0;
switch normType
    case 'max'
        Fk = TopoParam(2,1)*(Kf./max(Kf.*~A,[],'all'));
    case 'sum'
        Fk = TopoParam(2,1)*(Kf./sum(Kf.*~A,'all'));
end

Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;

Ff = (Fd+Fk).*~A;

P = Ff(indx);


for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;
    EdgesTriu(b(i)) = 1;
    uu = u(r);
    vv = v(r);
    x = A(uu,:);
    y = A(:,vv);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    K(uu,y) = K(uu,y) + 1;
    K(y,uu) = K(y,uu) + 1;
    K(vv,x) = K(vv,x) + 1;
    K(x,vv) = K(x,vv) + 1;

    switch mv2
        case 'powerlaw'
            Fk = K.^TopoParam(1,1);
        case 'exponential'       
            Fk = exp(TopoParam(1,1)*K);    
    end
    Fk(isinf(Fk)) = 0;
    
    switch normType
        case 'max'
            for j = 1:length(PDMsmv)
                Df = Df1(:,:,j);
                PDf(:,:,j) = PDMsParam(2,j)*(Df./max(Df.*~A,[],'all'));
            end
            Fd = sum(PDf,3);
            MaxNorm = max(Fk.*~A,[],'all');
            if MaxNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = TopoParam(2,1)*(Fk./MaxNorm).*~A;   
            end
        case 'sum'
            for j = 1:length(PDMsmv)
                Df = Df1(:,:,j);
                PDf(:,:,j) = PDMsParam(2,j)*(Df./sum(Df.*~A,'all'));
            end
            SumNorm = sum(Fk.*~A,'all')/2;
            if SumNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = TopoParam(2,1)*(Fk./SumNorm).*~A;   
            end
    end
        
    Ff = Fd + TopoTerm;

    Ff(A) = 0;

    P = Ff(indx);
    
end
b = indx(b);

function b = fcn_matching(A,K,PDMs,m,PDMsParam,TopoParam,modelvar,epsilon,normType)
% Convert input matrix A to a logical
A = A > 0;

% Add epsilon to matrix K, where K is the inital matching matrix

K = K + epsilon;

% Get the length of vector Fd. Fd is the distance matrix with the
% powerlaw/exponential decay already applied
n = length(A);

% Calculate the number of edges to add (non-zero elements in A divided by 2)
mseed = nnz(A)/2;

% Find the upper triangular indices of a ones matrix
[u,v] = find(triu(ones(n),1));

% Calculate indices for an n*n matrix
indx = (v - 1)*n + u;

% Initialize output vector b
b = zeros(m,1);

% Assign the indices of existing edges in b
b(1:mseed) = find(A(indx));

EdgesTriu = zeros(size(indx));
EdgesTriu(b(b~=0)) = 1;

PDMsmv = modelvar{1};
PDf = zeros(size(PDMs));
Df1 = zeros(size(PDMs));
for j = 1:length(PDMsmv)
    mv1 = PDMsmv{j};
    switch mv1
        case 'powerlaw'
            Df1(:,:,j) = PDMs(:,:,j).^PDMsParam(1,j);
        case 'exponential'       
            Df1(:,:,j) = exp(PDMsParam(1,j)*(PDMs(:,:,j))).*~eye(n);                 
    end
    Df = Df1(:,:,j);
    switch normType
    case 'max'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./max(Df.*~A,[],'all'));
    case 'sum'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./sum(Df.*~A,'all'));        
    end
end

Fd = sum(PDf,3);

mv2 = modelvar{2};

switch mv2
    case 'powerlaw'
        Fk = K.^TopoParam(1,1);
    case 'exponential'       
        Fk = exp(TopoParam(1,1)*K);    
end
Fk(isinf(Fk)) = 0;
switch normType
    case 'max'
        TopoTerm = TopoParam(2,1)*(Fk./max(Fk(indx(~EdgesTriu))));
    case 'sum'
        TopoTerm = TopoParam(2,1)*(Fk./sum(Fk(indx(~EdgesTriu))));
end

TopoTerm(isnan(TopoTerm)) = 0;TopoTerm(isinf(TopoTerm)) = 0;
    
Fp = (Fd + TopoTerm).*~A;

% Get probabilities P based on Ff and indx (extracts each edge only once)
P = Fp(indx);

% Loop to find additional edges
for ii = (mseed + 1):m
    
    
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    
    b(ii) = r;

    EdgesTriu(b(ii)) = 1;
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
            Fk = K.^TopoParam(1,1);
        case 'exponential'       
            Fk = exp(TopoParam(1,1)*K);    
    end
    Fk(isinf(Fk)) = 0;  
        
    switch normType
        case 'max'
            for j = 1:length(PDMsmv)
                Df = Df1(:,:,j);
                PDf(:,:,j) = PDMsParam(2,j)*(Df./max(Df.*~A,[],'all'));
            end
            Fd = sum(PDf,3);
            MaxNorm = max(Fk.*~A,[],'all');
            if MaxNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = TopoParam(2,1)*(Fk./MaxNorm).*~A;   
            end
        case 'sum'
            for j = 1:length(PDMsmv)
                Df = Df1(:,:,j);
                PDf(:,:,j) = PDMsParam(2,j)*(Df./sum(Df.*~A,'all'));
            end
            SumNorm = sum(Fk.*~A,'all')/2;
            if SumNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = TopoParam(2,1)*(Fk./SumNorm).*~A;   
            end
    end
       
    Fp = (Fd + TopoTerm).*~A;    
    % Update probabilities P based on Ff and indx
    P = Fp(indx);

end
b = indx(b);

function b = fcn_sptl(A,PDMs,m,PDMsParam,modelvar,normType)
n = length(A);
mseed = nnz(A)/2;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;

PDMsmv = modelvar{1};
PDf = zeros(size(PDMs));
Df1 = zeros(size(PDMs));
for j = 1:length(PDMsmv)
    mv1 = PDMsmv{j};
       
    switch mv1
        case 'powerlaw'
            Df1(:,:,j) = PDMs(:,:,j).^PDMsParam(1,j);
        case 'exponential'       
            Df1(:,:,j) = exp(PDMsParam(1,j)*(PDMs(:,:,j))).*~eye(n);                 
    end
    Df = Df1(:,:,j);
    switch normType
    case 'max'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./max(Df.*~A,[],'all'));
    case 'sum'
        PDf(:,:,j) = PDMsParam(2,j)*(Df./sum(Df.*~A,'all'));        
    end
end

Fd = sum(PDf,3);

P = Fd(indx).*~A(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;

    EdgesTriu(b(i)) = 1;
    
    for j = 1:length(PDMsmv)
       Df = Df(:,:,j);
       switch normType        
        case 'max'
            PDf(:,:,j) = PDMsParam(2,i)*(Df./max(Df.*~A,[],'all'));
        case 'sum'
            PDf(:,:,j) = PDMsParam(2,i)*(Df./sum(Df.*~A,'all'));      
       end
    end
    Fd = sum(PDf,3);
    P = Fd(indx);
    P(b(1:i)) = 0;
end
b = indx(b);