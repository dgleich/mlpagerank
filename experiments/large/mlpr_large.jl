function readSMAT(filename)
    (rows,header) = readdlm(filename;header=true)
    A = sparse(round(Int64,rows[:,1])+1, round(Int64,rows[:,2])+1, rows[:,3], 
               parse(Int64,header[1]), parseint(Int64,header[1]))
    return A
end

function test_tricycles()
    A = [0 1 0; 1 0 1; 1 1 0];    
    A = float64(A)
    A = sparse(A)
    
    return tricycles(A)
end


function tricycles(G::SparseMatrixCSC{Float64,Int64})

    Gt = G'
    n = G.n
    
    #tris = zeros(Int64,0,3)
    trilist = zeros(Int64,0)
    
    outcolptr = Gt.colptr
    outrowval = Gt.rowval
    incolptr = G.colptr
    inrowval = G.rowval
    
    neighs = zeros(Int32,n)
    
    for es=1:n
        for nzi in outcolptr[es]:(outcolptr[es+1] - 1)
            ed = outrowval[nzi]   
            
            # so we have an edge (es,ed)
            
            # index out-neighs of ed
            
            for nzi2 = outcolptr[ed]:(outcolptr[ed+1] - 1)
                neighs[outrowval[nzi2]] = 1
            end
            
            # look at in-edges to es
            for nzi3 = incolptr[es]:(incolptr[es+1] - 1)
                vt = inrowval[nzi3]
                if neighs[vt] > 0
                    push!(trilist, es)
                    push!(trilist, ed)
                    push!(trilist, vt)
                end
            end
            
            # clear out-neighs index
            for nzi2 = outcolptr[ed]:(outcolptr[ed+1] - 1)
                neighs[outrowval[nzi2]] = 0
            end
        end
    end
    
    tris = reshape(trilist, 3, int64(length(trilist)/3))
    return tris
end

function get_norms(n::Int64,tris::Array{Int64,2})
    ntris = size(tris,2)
    norms = zeros(Int64,ntris)
    
    # todo, could optimize this step, but it's pretty good
    M = sparse(vec(tris[2,:]),vec(tris[3,:]),1,n,n)
    
    for i=1:ntris
        norms[i] = M[tris[2,i],tris[3,i]]
    end
    
    return norms, M
end

function trimult(n::Int64,tris::Array{Int64,2},norms::Array{Int64,1},x::Array{Float64,1})
    y = zeros(Float64,n)
    
    for ti=1:size(tris,2)
        ci = tris[1,ti]
        cj = tris[2,ti]
        ck = tris[3,ti]
        
        y[ci] += x[cj]*x[ck]/norms[ti]
    end
    
    return y
end

function mlpr_cycles(n::Int64, A::SparseMatrixCSC{Float64,Int64}, 
    tris::Array{Int64,2}, alpha::Float64;
    tol::Float64=1e-16, shift::Float64=0., maxiter::Int64=1000,
    beta::Float64=0.5)
    verbose=true
    
    # compute the data we need
    (colnorms,M) = get_norms(n,tris)
    degs = sum(A,2)
    
    v = ones(n)/n
    x = copy(v)
    
    hist=zeros(Float64,maxiter)
    niter = 0
    
    for iter=1:maxiter
        # we are going to use the following three-way tensor
        # 0.5*Markov M + 0.5*triangles
        # 
        y = trimult(n,tris,colnorms,x);
        sumy = sum(y)
        y = y + (1-sumy)*v
        
        
        z = A'*(vec(x)./vec(degs))
        sumz = sum(z)
        z = z + (1-sumz)*v

        #y = y + (1-sumy)*z
        #y = z + 0.

        y = beta*y + (1-beta)*z       
 
        y = alpha*y + (1-alpha)*v;
    
        hist[iter] = norm(y-x,1);
                
        if verbose
            @printf("iter %5i  resid %8.2e  sumy %8.2e\n", iter, hist[iter], sumy)
        end
    
        x = (1-shift)*y + shift*x;
        niter = iter
        
        if hist[iter] < tol
            break
        end

    end
    
    hist = hist[1:niter]
    
    return x, hist
end
    
    

function mlpr_cycles(G::SparseMatrixCSC{Float64,Int64},alpha::Float64, beta::Float64)
    n = G.n
    
    # compute the data we need
    tris = tricycles(G)

    return mlpr_cycles(n,A,tris,alpha;beta=beta)    
end
    
A = readSMAT("wiki-Talk.smat")
@printf("Done reading...\n");
(x,hist) = mlpr_cycles(A,0.99,0.5)

alphas=[0.5 0.85 0.99]
betas=[0.5 0.85 0.99]
niter = zeros(length(alphas),length(betas))

for ai=eachindex(alphas)
    for bi=eachindex(betas)
        alpha = alphas[ai]
        beta = betas[bi]
        (x,hist) = mlpr_cycles(A,alpha,beta)
        niter[ai,bi] = length(hist)
    end
end

@printf("n nodes: %i\n", size(A,1))
@printf("nnz in tricycles: %i\n",maximum(size(tricycles(A))))


for ai=eachindex(alphas)
    @printf("%.2f & %s \\\\ \n", alphas[ai],
        join([@sprintf("%i",x) for x in niter[ai,:]]," & "))
end
