function z=LaguerreGaussianE(params,r,varargin);

if nargin>=3
    if isstr(varargin{1})
        if strcmp(varargin{1},'cart')          % use cartesian coordinates
            defaultcoord2=1;
            cartesianflag=1;
        else                                    % use polar coordinates with the default theta
            defaultcoord2=1;
            cartesianflag=0;
        end
    else                                        % use polar coordinates with the specified theta
        defaultcoord2=0;
        cartesianflag=0;
        theta=varargin{1};
    end
else                                            % use polar coordinates with the default theta              
    defaultcoord2=1;
    cartesianflag=0;
end

if nargin>=4
    defautlcoord2=0;                    
    if isstr(varargin{2}) & strcmp(varargin{2},'cart')
        cartesianflag=1;                         % use cartesian coordinates with specified y
    else                                         % use polar coordinates with the specified theta
        cartesianflag=0;
    end
end


if cartesianflag                                                    % cartesian (x,y) domain supplied
    x=r;
    if min(size(x))==1                                              % map is 2->1 on a cartesian domain
        if size(x,1)<size(x,2), x=x'; end                           % make x and y columnar
        if defaultcoord2
            y=zeros(size(x));
        else 
            y=theta;
            if size(y,1)<size(y,2), y=y'; end
        end 
    end
    if min(size(x)) > 1                                             % map is 2->2 on a cartesian domain
        if defaultcoord2 
            y=transpose(x); 
        else
            y=theta;
        end
        z=zeros(size(x,1),size(x,2),size(params,1));                % need this since zeros(size(y),10) gives a 2D matrix even if y is 2D!  (Matlab feature.)
    else
        z=zeros(size(x),size(params,1)); 
    end
    [theta,r]=cart2pol(x,y);                                        % convert to polar coords for calculation
else                                                                % polar (r,theta) domain supplied
    if min(size(r))==1                                              % map is 2->1 on a polar domain
        if size(r,1)<size(r,2), r=r.'; end                          % make r columnar
        if defaultcoord2                                            % make theta columnar
            theta=zeros(size(r));                                   % default 1D theta is zero
        else 
            if size(theta,1)<size(theta,2), theta=theta.'; end
        end
    else                                                            % otherwise assume r and theta are already in meshgrid format
        z=zeros(size(r,1),size(r,2),size(params,1));                % need this since zeros(size(r),10) gives a 2D matrix even if y is 2D!  (Matlab feature.)
    end
end


p=params(:,1);
m=params(:,2);
signm=sign(m);
m=abs(m);
q=params(:,3);
if size(params,2)>=4
    lambda=params(:,4);
else
    lambda=1064e-9;
end
if size(params,2)>=5
    a=params(:,5);
else
    a=ones(size(q));
end

w=w_(q,lambda);

if min(size(r))>=2
    %{
    for u=1:size(params,1)
            z(:,:,u) = a(u)...
                .* sqrt(2*factorial(p(u))/(1+(m(u)==0))/pi/(factorial( m(u)+p(u) )))/w(u)...
                .* (sqrt(2)*r/w(u)).^m(u) .*exp(i*signm(u)*m(u).*theta).* LaguerrePoly([p(u),m(u)],2*r.^2/w(u).^2)...
                .* exp( -i*2*pi/lambda(u)*r.^2/2/q(u));
    end
      %}
    for u=1:size(params,1)
        z(:,:,u)=...    
                (1/w(u)).*sqrt(2*factorial(p(u))/(pi.*factorial(m(u)+p(u))))...
            .*exp(-r.^2/w(u).^2).*(sqrt(2)*r/w(u)).^m(u).* LaguerrePoly([p(u),m(u)],2*r.^2/w(u).^2)...
            .*exp(+i*m(u)*theta);      
    end
  
else
    for u=1:size(params,1)
            z(:,u) = a(u)...
                .* sqrt(2*factorial(p(u))/(1+(m(u)==0))/pi/(factorial( m(u)+p(u) )))/w(u)...
                .* (sqrt(2)*r/w(u)).^m(u) .* exp(i*signm(u)*m(u).*theta).* LaguerrePoly([p(u),m(u)],2*r.^2/w(u).^2)...
                .* exp( -i*2*pi/lambda(u)*r.^2/2/q(u));
    end
end