% Vector form of factorial.
% For some reason the function "factorial" can't act on vectors.
% 
function y=fact(x);

y=[];
for s=1:length(x)
    y(s)=factorial(x(s));
end