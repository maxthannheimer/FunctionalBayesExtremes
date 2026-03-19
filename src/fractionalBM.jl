""" Function to generate a fractional Brownian motion samples """


using FFTW
using Distributions
using LinearAlgebra


""" Standard way to generate fBm samples, should work for small Grids. For later comparison and coarse grid simulations. """
function fBm(;param::Parameter,grid::Grid,num_sim::Int)::Vector{Vector{Float64}}
    c_sqrt=sqrt(param.c)
    # Generate the covariance matrix
    N = grid.gridsize^2   
    cov_matrix = zeros(Float64, N, N)
    tx = grid.coord_fine  #grid points in x and y
    for i in 1:N
        for j in 1:N
            cov_matrix[i, j] =  (norm(grid.coord_fine[i,:] )^param.β + norm(grid.coord_fine[j,:] )^param.β - norm(grid.coord_fine[i,:] - grid.coord_fine[j,:])^param.β).+5*1e-6 
        end
    end

    # Generate a sample from the multivariate normal distribution
    mean_vector = zeros(Float64, N)
    #return cov_matrix
    res1=[Vector{Float64}(undef,N) for i in 1:num_sim]
    for i in 1:num_sim
        res1[i] = rand(MvNormal(mean_vector, cov_matrix))
    end
    return c_sqrt.*res1
end

""" Analog to fft2 function from matlab, but for 2D arrays. """
function fft2d(arr)
    return fft(fft(arr, 1), 2)
end     

function fft2(A)
    # Apply 1D FFT along each dimension
    fft_rows = fft(A, 1)
    fft_cols = fft(fft_rows, 2)

    # Return the result
    return fft_cols
end


""" modified cov function """
function rho(x,y,R,a)
    #embedding of covariance function on a larger [0,R] × [0,R] Grid
    if a<1.5 #alpha=2 Hurst param TODO: check if this is correct, literature says a<=1.5, we used a<=1.1
        beta=0
        c_2=a/2
        c_0=1-a/2
    else
        #TODO: R=2 in literature table, but no R=2 in code, check whats right 
        beta=a*(2-a)/(3*R*(R^2-1))
        c_2=(a-beta*(R-1)^2*(R+2))/2
        c_0=beta*(R-1)^3+1-c_2
    end
    #create cont isotropic cov function
    r=sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)
    if r<=1
        out=c_0-r^a+c_2*r^2
   #elseif R==1 #caution?
    #    out=0
    elseif r<=R
        out=beta*(R-r)^3/r
    else
        out=0
    end
    return (out,c_0,c_2)
end

""" Function to generate a fractional Brownian motion samples using circulant embedding method """

function FBM_simu_fast(;param::Parameter,grid::Grid,num_sim::Int)::Vector{Vector{Float64}}
    gridsize=grid.gridsize
    c_sqrt=sqrt(param.c)
    H=param.β/2 #Hurst Param
    if param.β<1.5 #TODO: check if this is correct, literature says a<=1.5, we used a<=1.1
        R=1
    else
        R=2 #expanded gridsize, region of interest is [0,1]×[0.1]
    end
    n=m=R*gridsize #size of grid is m ∗ n, cov mat size is n^2*m^2
    tx=(1:n)/n*R; ty=(1:m)/m*R #grid points in x and y 
    Rows=zeros(m,n)
    for i in 1:n 
        for j in 1:m
            Rows[j,i]=rho([tx[i],ty[j]],[tx[1],ty[1]],R,2*H)[1]
        end
    end
    BlkCirc_row=[Rows  Rows[:,end-1:-1:2] ;  Rows[end-1:-1:2,:]  Rows[end-1:-1:2,end-1:-1:2 ]  ]
    #calculate eigenvalues via fft
    eig_vals=real(fft2(BlkCirc_row)/(4*(m-1)*(n-1)))
    #optional:
    #set small values to zero:
    #eps = 10^-8
    #eig_vals=eig_vals.*(abs.(eig_vals) .> eps)
    #eig_vals=eig_vals.*(eig_vals .>= 0)
    eig_vals=sqrt.(eig_vals)
   
    res1=[Vector{Float64}(undef,gridsize*gridsize) for i in 1:num_sim]
    #one can get two times as many simulations for free by using the imaginary and real part 
    #of the complex gaussian, but they are dependend
    
    #res2=[Matrix{Float64}(undef,gridsize,gridsize) for i in 1:num_sim]
    
    for trial in 1:num_sim
        #generate field with covariance given by block circulant matrix
        Z= randn(2*(m-1),2*(n-1)) + im* randn(2*(m-1),2*(n-1)) 
        #fft to calculate diag*Z
        F=fft2(eig_vals.*Z)
        #extract subblock with desired cov variance
        F=F[1:gridsize,1:gridsize]
        (out,c_0,c_2)=rho([0,0],[0,0],R,2*H)
        #optional two (dependend) real fields
        field1=real(F)
        #field2=imag(F)
        #set field zero at origin
        field1=field1.-field1[1,1]
        #field2=field2.-field2[1,1]
        #make correction for embedding with a term c_2*r^2
        X_grid_comp = [i for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
        Y_grid_comp = [j for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
    
        res1[trial]=vec((field1 + (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2))')
        #res2[trial]=field2 +  (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2)
    end

    #(c_sqrt.*res1,c_sqrt.*res2)
    c_sqrt.*res1
end



#param=Parameter(α=1.0, β=2.0, c=3.0)
#grid=Grid()

#num_sim=100
#FBM_res=FBM_simu_fast(param=param,grid=grid,num_sim=num_sim)
#FBM_covmat_res=fBm(param=param,grid=grid,num_sim=num_sim)
