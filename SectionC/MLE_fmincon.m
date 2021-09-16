clear all 
n=2; 
thetatrue=[1 1 2 2 0.2]; 
mu = [0 0]; %mean and variance covariance matrix
sd = [1 thetatrue(2*n+1); thetatrue(2*n+1) 1];
data = randi([0 8],100,4); % make integer array of 100 observations of 4 variables
X=data(:,1:n); % x values are first n columns of array
Y=data(:,n+1:size(data,2));  % y values are rest of columns
B1(:,2)=-X(:,1); 
B1(:,1)=-1;   
B2(:,2)=-X(:,2);
B2(:,1)=-1;
C1=(all(bsxfun(@eq,Y,[0 0]),2)); % change to 0s and 1s using eq function and return logical true if all values in second dimension are nonzero/true
C2=1-C1;   % return 0 or 1, depending on previous result
cdf=@(x) mvncdf( [B1*[x(1);x(3)],  B2*[x(2);x(4)] ] ,mu,[1 x(5); x(5) 1]);             
options=optimset('Algorithm',...
'interior-point','Display','iter','MaxIter',10000,'TolX',10^-30,'TolFun',10^-30);
theta0=thetatrue; % set initial theta randomly
[theta,fval,exitflag,output]=...
        fmincon(@(x) log_lik(x,cdf,C1,C2),theta0,[],[],[],[],[-Inf; -Inf; -Inf; -Inf; -1], ...
                 [+Inf; +Inf; +Inf; +Inf; 1],[],options);