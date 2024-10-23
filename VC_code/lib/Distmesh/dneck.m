function d=dneck(p,x1,x2,h1,h2)

d=-min(min(min((h1-h2)/2*cos((p(:,1)-x1)/(x2-x1)*pi)+(h1+h2)/2+p(:,2),(h1-h2)/2*cos((p(:,1)-x1)/(x2-x1)*pi)+(h1+h2)/2-p(:,2)),-x1+p(:,1)),x2-p(:,1));
