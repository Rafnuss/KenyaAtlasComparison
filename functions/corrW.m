function corr = corrW(x,y,w)

% deal with nan
id = isnan(x)|isnan(y)|isnan(w);
id=~id(:);
x=x(id);
y=y(id);
w=w(id);
w = w./sum(w);

mx = w'*x;
my = w'*y;

sx = w'*(x-mx).^2;
sy = w'*(y-my).^2;

sxy = w'*((x-mx).*(y-my));

corr = sxy / sqrt( sx * sy);

end