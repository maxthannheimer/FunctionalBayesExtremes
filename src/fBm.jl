""" Function to generate a fractional Brownian motion samples """


using FFTW


""" AI generated code, may contain errors. Maybe not the most efficient way to generate fBm samples, but it should work for small N. For later comparion """
function fBm(H::Float64, N::Int)
    # Generate the covariance matrix
    cov_matrix = zeros(Float64, N, N)
    for i in 1:N
        for j in 1:N
            cov_matrix[i, j] = 0.5 * (abs(i)^H + abs(j)^H - abs(i - j)^H)
        end
    end

    # Generate a sample from the multivariate normal distribution
    mean_vector = zeros(Float64, N)
    return rand(MvNormal(mean_vector, cov_matrix))
end

""" Analog to fft2 function from matlab, but for 2D arrays. """
function fft2d(arr::Matrix{Float64})::Matrix{ComplexF64}
    return fft(fft(arr, 1), 2)
end     

function fft2(A::Matrix{Float64})::Matrix{ComplexF64}
    # Apply 1D FFT along each dimension
    fft_rows = fft(A, 1)
    fft_cols = fft(fft_rows, 2)

    # Return the result
    return fft_cols
end


A=rand(4,4)
fft2d(A)-fft2(A)