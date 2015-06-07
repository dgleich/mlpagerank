module Puma

export readSMAT

function readSMAT{Ti,Tv}(filename;IntType::Type{Ti}=Int64,ValType::Type{Tv}=Float64)
    (rows,header) = readdlm(filename;header=true)
    
    A = sparse(round(IntType,rows[:,1])+1, round(IntType,rows[:,2])+1, 
                rows[:,3], 
               parseint(IntType,header[1]), parseint(IntType,header[1]))
    return A
end

# need components

# need normout
function normout!(A::SparseMatrixCSC)
    rowsums = sum(A,2)
    nzval = nonzeros(A)
    nzval ./
end