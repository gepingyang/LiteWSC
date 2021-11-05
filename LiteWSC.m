%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the code for the LiteWWSC algorithm, which is proposed in   %
% the following paper:                                                %
%LiteWSSC:a Lightweight framework for Web-Scale Spectral Clustering   %                                             %
%                                                                     %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [label_orig, label] = LiteWWSC(dataname, k, s, p, r)

if nargin < 5
    r = 5; % The number of nearest prototypes.
end
if nargin < 4
    p = 500; % The number of prototypes.
end
if nargin < 3
    s = 0; % tthe number of data points that are used to generate prototypes and prototype graph.
            % s = 0 denotes that  all data of this btach is used to generate prototypes and
            % prototype graph.
end 
if nargin < 2 
    disp('you have to input the number of target clustering');
    return;
end
if nargin < 1 | size(dataname) < 1
    disp('you have to input data');
    return;
end
flg = 1;% flg = 1 denotes the first time loading data, we sholud contrcut prototypes graph, otherwise, we just need to assign data points to
        % their nearest prototpyes.
for dataname_batch = dataname
     if flg == 1
         flg = 0;
        %loading adta from hard disk 
        dataname_batch = char(dataname_batch);
        load( dataname_batch,'fea','gnd');
       N = size(fea, 1);
        fea = full(fea); 
        label_orig = gnd;%sotre the original label of data
 %%%%%%%%%%%%%%%%%find prototypes%%%%%%%%%%%%%%%%%%%%       
        if  s == 0
            [label1, prototypes] = litekmeans(fea, p ,'MaxIter', 3,'Replicates',1);
        else
            indSmp = randperm(N);
            [label1, prototypes] = litekmeans(fea(indSmp(1:s),:), p ,'MaxIter', 3,'Replicates',1);
        end
       clear label1

%%%%%%%prototype graph construction%%%%%%%%%%%%%%%%%%%%%%
       D = EuDist2(fea,prototypes,0);
       [dump1, idx1] = min(D, [], 2);
       D = D(indSmp(1:s),:);
       sigma = mean(mean(D));
       dump = zeros(s,r);
       idx = dump;
       dump(:, 1) = dump1(indSmp(1:s));
       clear dump1
       idx(:,1) = idx1(indSmp(1:s));

       for ii = 2:r
           temp = (idx(:,ii - 1)-1)*s+[1:s]';
           D(temp) = 1e100; 
           [dump(:,ii),idx(:,ii)] = min(D,[],2);
       end
       clear D
       dump = exp(-dump/2/sigma^2);
       Gsdx = dump;
       Gidx = repmat([1:s]',1,r);
       Gjdx = idx;
       clear dump 
       Z=sparse(Gidx(:),Gjdx(:),Gsdx(:),s,p);
       A = Z'*Z; %prototypes affinity matrix
       D_p = 1./(sqrt(sum(A,2))+10^(-6));%  prototypes degree matrix
       L = repmat(D_p,1,p).*A.*repmat(D_p', p, 1); %prototypes Laplacian matrix
       clear A
%%%%%%% %%%eigendedecomposition%%%%%%%%%%%%%%%%
       if issparse(L)
            L = full(L);
       end
        L = (L+L')/2;
        [U, eigvalue] = eig(L);
        eigvalue = diag(eigvalue);
        [dump, index] = sort(-eigvalue);
        clear dump
        U = U(:, index);
        U = U(:,2: k+1) ;
        U = U ./repmat(sqrt(sum(U.^2,2)),1,k);
 %%%%%%%%%%%%%%%%peroforming $k$-means%%%%%%%%%%%%%%%%
        [labels]=litekmeans(U,k,'MaxIter',100,'Replicates',10);
        %assign data points to their nearest prototypes.
        label = labels(idx1);
        clear  U  fea gnd  L D_a A
            
     else
         %load data in batch and assign data  to their nearest prototypes.
         dataname_batch = char(dataname_batch);
         load( dataname_batch,'fea', 'gnd');
         fea = full(fea);   
         label_orig = [label_orig;gnd];
         D = EuDist2(fea, prototypes, 0);
         [dump, idx2] = min(D, [], 2);
         label_batch = labels(idx2);
         label = [label;label_batch];
         clear  dump D fea gnd dump idx2 label_batch
     end
          
end
end


