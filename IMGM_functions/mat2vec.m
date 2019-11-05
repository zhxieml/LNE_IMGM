function p = mat2vec(X)%column-wise vectorization
    %row-wise vectorization
     X=X';
    p=X(:);
end