
function llr = mcmc_eq2(x, h, lambda2e, mc_iter, nvar)

%Generate init sample bit sequence randomly
L = length(h);
N = length(lambda2e);

p1 = 1./(1+exp(lambda2e));
p0 = 1-p1;

candidate = (rand(1, N)>p0);

sa = 1-2*candidate;
sa = [zeros(1,L-1), sa, zeros(1,L-1)];

%Run Gibbs sampler to generate bit sequence table
bitsTable = zeros(mc_iter, N);

for ii=1:mc_iter
    for n=1:N
        tmp0=0;
        tmp1=0;
        
        for k=n:n+L-1
          sa0 = sa(k+L-1:-1:k);
          sa0(k-n+1) = 1;
          sa1 = sa(k+L-1:-1:k);
          sa1(k-n+1) = -1;
          tmp0 = tmp0+(x(k)-sum(h.*sa0))^2;  
          tmp1 = tmp1+(x(k)-sum(h.*sa1))^2;  
        end
        gama0 = exp(-1/(2*nvar)*tmp0+lambda2e(n)/2);
        gama1 = exp(-1/(2*nvar)*tmp1-lambda2e(n)/2);
        pp0 = gama0/(gama0+gama1);
        newbit = (rand(1)>pp0);
        candidate(n) = newbit;
        sa(n+L-1) = 1-2*newbit;
    end
    bitsTable(ii,:) = candidate;
end

% compute LLR
bitTableLen = mc_iter;
for n=1:N
    max0 = -1.0e7;
    max1 = -1.0e7;
    
    if n<L %first L-1 bits
        sa = zeros(1,L-n); 
        %Get parital table
        In1_n2 = bitsTable(:,1:n+L-1); 
        %flip bit n to get extending table
        tmpTab = In1_n2;
        tmpTab(:,n) = ~In1_n2(:,n); 
        In1_n2 = [In1_n2; tmpTab];
        
        %remove repeated bit sequence
        In1_n2_nrep = zeros(1, size(In1_n2, 2));
        In1_n2_nrep(1,:) = In1_n2(1,:);
        nrepBitTableLen = 1;
        for jj=2:2*bitTableLen
            repeatFlag = 0;
            for kk=1:nrepBitTableLen
                if sum(abs(In1_n2(jj,:)-In1_n2_nrep(kk,:)))==0
                    repeatFlag = 1;
                    break;
                end
            end
            if repeatFlag==0
                nrepBitTableLen = nrepBitTableLen+1;
                In1_n2_nrep(nrepBitTableLen,:) = In1_n2(jj,:);
            end    
        end
        
        %Calculate max log-probability to be 0 or 1
        for jj=1:nrepBitTableLen  
            sa1(1:L-n) = sa;
            sa1(L-n+1:2*L-1) = 1-2*In1_n2_nrep(jj,:);
            tmp = 0;
            for k=1:L
                sa2 = sa1(k+L-1:-1:k);
                tmp = tmp+(x(n+k-1)-sum(h.*sa2))^2;
            end
            
            tmp = -tmp/(2*nvar);
            for k=1:n+L-1
                if In1_n2_nrep(jj,k)==0
                    tmp = tmp + lambda2e(k);
                end
            end
            
             if In1_n2_nrep(jj,n)==0
                 if tmp>max0
                     max0=tmp;
                 end
             else
                 if tmp>max1
                     max1=tmp;
                 end
             end
        end
    elseif n<N-L % In the middle bits
        In1_n2 = bitsTable(:,n-L+1:n+L-1); 
        %flip bit n to get extending table
        tmpTab = In1_n2;
        tmpTab(:,L) = ~In1_n2(:,L); 
        In1_n2 = [In1_n2; tmpTab];
        
        %remove repeated bit sequence
        In1_n2_nrep = zeros(1, size(In1_n2, 2));
        In1_n2_nrep(1,:) = In1_n2(1,:);
        nrepBitTableLen = 1;
        
        for jj=2:2*bitTableLen
            repeatFlag = 0;
            for kk=1:nrepBitTableLen
                if sum(abs(In1_n2(jj,:)-In1_n2_nrep(kk,:)))==0
                    repeatFlag = 1;
                    break;
                end
            end
            if repeatFlag==0
                nrepBitTableLen = nrepBitTableLen+1;
                In1_n2_nrep(nrepBitTableLen,:) = In1_n2(jj,:);
            end    
        end
        
        %Calculate max log-probability to be 0 or 1
        for jj=1:nrepBitTableLen  
            sa1(1:2*L-1) = 1-2*In1_n2_nrep(jj,:);
            tmp = 0;
            for k=1:L
                sa2 = sa1(k+L-1:-1:k);
                tmp = tmp+(x(n+k-1)-sum(h.*sa2))^2;
            end
            
            tmp = -tmp/(2*nvar);
            for k=1:2*L-1
                if In1_n2_nrep(jj,k)==0
                    tmp = tmp + lambda2e(k);
                end
            end
            
            if In1_n2_nrep(jj,L)==0
                 if tmp>max0
                     max0=tmp;
                 end
             else
                 if tmp>max1
                     max1=tmp;
                 end
             end
        end
    else %tail bits
        In1_n2 = bitsTable(:,n-L+1:N); 
        %flip bit n to get extending table
        tmpTab = In1_n2;
        tmpTab(:,L) = ~In1_n2(:,L); 
        In1_n2 = [In1_n2; tmpTab];
        
        %remove repeated bit sequence
        In1_n2_nrep = zeros(1, size(In1_n2, 2));
        In1_n2_nrep(1,:) = In1_n2(1,:);
        nrepBitTableLen = 1;
        
        for jj=2:2*bitTableLen
            repeatFlag = 0;
            for kk=1:nrepBitTableLen
                if sum(abs(In1_n2(jj,:)-In1_n2_nrep(kk,:)))==0
                    repeatFlag = 1;
                    break;
                end
            end
            if repeatFlag==0
                nrepBitTableLen = nrepBitTableLen+1;
                In1_n2_nrep(nrepBitTableLen,:) = In1_n2(jj,:);
            end    
        end
        
        %Calculate max log-probability to be 0 or 1
        for jj=1:nrepBitTableLen  
            sa1(1:N-n+L) = 1-2*In1_n2_nrep(jj,:);
            sa1(N-n+L+1:2*L-1)=0;
            tmp = 0;
            for k=1:L
                sa2 = sa1(k+L-1:-1:k);
                tmp = tmp+(x(n+k-1)-sum(h.*sa2))^2;
            end
            
            tmp = -tmp/(2*nvar);
            for k=1:N-n+L
                if In1_n2_nrep(jj,k)==0
                    tmp = tmp + lambda2e(k);
                end
            end
            
            if In1_n2_nrep(jj,L)==0
                 if tmp>max0
                     max0=tmp;
                 end
             else
                 if tmp>max1
                     max1=tmp;
                 end
             end
        end 
    end
        
    llr(n) = max0-max1-lambda2e(n); %max-log 
end