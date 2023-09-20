function [B,b,t] = matching_gen_model_mult_old(A,PD,m,modelvar,PDexpo,gam,epsilon)
%           Run generative model code for the multiplicative model
% 
%   Generates synthetic networks using the models described in the study by
%   Betzel et al (2016) in Neuroimage.
%
%   Inputs:
%           A,          binary network of seed connections
%           PD,         Euclidean distance/fiber length/node similarity
%                       matrix. Multiple can be input either as a cell, 
%                       where each cell contains a different matrix or as a
%                       3D matrix (n*n*nPD, where n is the number of nodes
%                       and nPD is the number of PD matrices).
%           m,          number of connections that should be present in
%                       final synthetic network
%           modelvar,   specifies whether the generative rules are based on
%                       power-law or exponential relationship
%                       ({'powerlaw'}|{'exponential})
%           PDexpo,     the parameter controlling the values in PD. If
%                       there are multipe PD matrices, PDexpo should be a
%                       vector where each index gives the parameter for the
%                       corresponding PD matrix
%           gam,        the parameter controlling topology
%           epsilon,    the baseline probability of forming a particular
%                       connection (should be a very small number
%                       {default = 1e-6}).
%
%   Output:
%           B,          an adjacency matrix
%           b,          a vector giving the index of each edge in B. Note
%                       that the ordering of b shows which edges formed
%                       first (e.g., b(1) was the first edge to form, b(2)
%                       the second etc etc).
%           t,          the time in seconds it took do do each iteration
%
%       How to convert b to B:
%       n = length(A); B = zeros(n); B(b(:,i)) = 1; B = B + B'; 
%
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%   Oldham et al (2022) Science Advances 10.1126/sciadv.abm6127
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015
%   Edited by Stuart Oldham, Monash University 2021, MCRI 2023

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