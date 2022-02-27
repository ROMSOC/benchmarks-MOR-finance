function [Y,I]=maxk(X,k)
    [Xs, Is]=sort(X);
    Y = Xs(end:-1:end-k+1);
    I = Is(end:-1:end-k+1);
end