%% Here we look at a graph based on what Nati Linial said
%
% We connect each node to anyone in it's "tube" area. So what happens if
% that graph is empty?

%% There is no 3x3 tensor that's fully stochastic that makes this work
n = 3; P = zeros(n,n,n); for i=1:n, P(i,i,:) = i; P(i,:,i) = i; P(:,i,i) = i; end
P = P == 0;
R = reshape(P,n,n^2);

%% There is no 4x4 tensor that's fully stochastic that makes this work
n = 4; P = zeros(n,n,n); for i=1:n, P(i,i,:) = i; P(i,:,i) = i; P(:,i,i) = i; end
P = P == 0;
R = reshape(P,n,n^2);

%% There is no 3x3 tensor that's fully stochastic that makes this work
n = 5; P = zeros(n,n,n); for i=1:n, P(i,i,:) = i; P(i,:,i) = i; P(:,i,i) = i; end
P = P == 0;
R = reshape(P,n,n^2);

%%
n = 6; P = zeros(n,n,n); for i=1:n, P(i,i,:) = i; P(i,:,i) = i; P(:,i,i) = i; end
P = P == 0;
R = reshape(P,n,n^2);