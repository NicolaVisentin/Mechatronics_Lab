function B=dfdu(u,data)

% extract data

C=data;

% compute dfdu

B=zeros(2,1);
B(2,1)=C*cos(u);

end