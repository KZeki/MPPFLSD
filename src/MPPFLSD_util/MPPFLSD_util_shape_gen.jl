function rms(u)
    return ( norm(u)/sqrt(length(u)) )
end

# function rms(u; n=2)
#     return ( sqrt(mean(u.^n)) )
# end

function diffmat(N::Int64, domain::AbstractArray, order::Int64)
    # swaping a and b so that, when the chebyshev points are calculated,
    # the domain is ordered from
    # small/left ((theta=pi)|>cos==-1) to large/right ((theta=0)|>cos==1)
    a = domain[2]
    b = domain[1]
    if N == 1
        return 0; # Trivial case for N=1
    elseif N <= 0
        return [;;]
    end
    
    # Compute Chebyshev-Gauss-Lobatto nodes in [a,b] 
    # no need to account for the constant part
    # since this gets cancled out in the D matrix anyways
    x_chebyshev =  (b - a) / 2 * cospi.((0:N-1) / (N-1)); # using cospi for accuracy

    # Compute first-order differentiation matrix for [-1,1] domain
    c_vec = [2; ones(N-2,1); 2] .* (-1).^(0:N-1);
    X = repeat(x_chebyshev,inner=(1,N));
    dX = X - X'; # Matrix of differences
    
    D_chebyshev = (c_vec * (1 ./ c_vec)')./(dX + I(N)); # Off-diagonal entries
    D_chebyshev = D_chebyshev - Diagonal(vec(sum(D_chebyshev,dims=2))); # Diagonal entries

    # Compute order differentiation matrix
    return (D_chebyshev)^order;
end
diffmat(N::Int64, order::Int64, domain::AbstractArray) = diffmat(N, domain, order)
diffmat(N::Int64) = diffmat(N, [-1, 1], 1)
diffmat(N::Int64, order::Int64) = diffmat(N, [-1, 1], order)
diffmat(N::Int64, domain::AbstractArray) = diffmat(N, domain, 1)

function chebpts(N,domain::AbstractArray=[-1,1])::Vector{Float64}
    # domain is defined as [a,b], but assigning [b,a]
    # so that the points are returned/scaled from smallest to largest
    # could also have used reverse() when returning,
    # but that is slower, even though only very slightly
    a = domain[2];
    b = domain[1];
    if N == 1 # Cheb point is a single point in the middle of the domain
        return (a + b) / 2
    elseif N <= 0
        return [;;]
    end
    x_chebpts = cospi.((0:N-1) / (N-1)); # using cospi for higher accuracy
    
    # Map to [a, b]
    x_chebpts_scaled = (b-a) / 2 * x_chebpts .+ (a + b) / 2; 
    
    return x_chebpts_scaled
end
        
function clencurt(N::Int64)::Vector{Float64} # clencurt
    #
    # Computes the integration weigths for pseudo-chebychev on domain [-1 1]
    #
    # INPUTS:
    # N  : the number of points 
    #
    # OUTPUT:
    # IW : vector of the integration weigths.
    
    nW=0:1:N-1;
    jW=0:1:N-1;
    
    bW=ones(1,N); 
    bW[1]=0.5; 
    bW[N]=0.5;
    cW=2*bW;
    bW=bW/(N-1);
    
    S=cos.(transpose(collect(nW[3:N])'.*collect(jW))*(pi/(N-1)));
    IW=bW.*( (2 .+(cW[3:N].*((1 .+(-1).^nW[3:N])./(1 .-nW[3:N].^2)))'*S) );
    
    return vec(IW)
end

function introw(N::Int64,domain::AbstractArray=[-1,1])::Vector{Float64}
    L = domain[2] - domain[1]
    w=L/2*clencurt(N);
    return w
end