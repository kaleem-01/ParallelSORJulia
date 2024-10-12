# using MPI
# using LinearAlgebra

# const N = 100           # Grid size (N x N)
# const MAX_ITER = 1000   # Maximum number of iterations
# const TOL = 1e-6        # Tolerance for convergence
# const OMEGA = 1.25      # Relaxation factor

# # Function to perform one SOR iteration on the local grid
# function sor_iteration(local_u, Nlocal, rank, nranks)
#     Nlocal_i, Nlocal_j = size(local_u)

#     for i in 2:Nlocal_i-1
#         for j in 2:Nlocal_j-1
#             local_u[i, j] = (1 - OMEGA) * local_u[i, j] + 
#                             OMEGA * 0.25 * (local_u[i-1, j] + local_u[i+1, j] + 
#                                             local_u[i, j-1] + local_u[i, j+1])
#         end
#     end
# end

# # Function to exchange boundaries between neighboring processes
# function exchange_boundaries(local_u, rank, size, Nlocal_i, Nlocal_j)
#     # Define neighboring ranks (assuming 1D decomposition for simplicity)
#     above = rank > 0 ? rank - 1 : MPI.PROC_NULL
#     below = rank < size - 1 ? rank + 1 : MPI.PROC_NULL

#     # Create send/receive buffers
#     send_top = local_u[2, :]  # Second row (since the first row is the boundary)
#     send_bottom = local_u[Nlocal_i-1, :]  # Second-to-last row (boundary row)
#     recv_top = similar(send_top)
#     recv_bottom = similar(send_bottom)

#     # Exchange boundaries with neighbors (above/below)
#     MPI.Sendrecv!(send_top, below, recv_bottom, above, MPI.COMM_WORLD)
#     MPI.Sendrecv!(send_bottom, above, recv_top, below, MPI.COMM_WORLD)

#     # Update local boundary rows after communication
#     if rank > 0
#         local_u[1, :] = recv_top  # Update top boundary
#     end
#     if rank < nranks - 1
#         local_u[Nlocal_i, :] = recv_bottom  # Update bottom boundary
#     end
# end

# # Main function to solve the Laplace equation
# function main()
#     MPI.Init()

#     rank = MPI.Comm_rank(MPI.COMM_WORLD)
#     nranks = MPI.Comm_size(MPI.COMM_WORLD)

#     # Grid decomposition: Each process gets a slice of the grid (1D decomposition)
#     Nlocal = N ÷ nranks
#     Nlocal_i = Nlocal + 2  # Adding 2 for the boundary rows
#     Nlocal_j = N  # Full columns, no decomposition in the x-direction

#     # Initialize the local grid
#     local_u = zeros(Float64, Nlocal_i, Nlocal_j)  # Solution array (including boundaries)

#     # Set boundary conditions (e.g., u = 1 on the left and right boundaries)
#     if rank == 0
#         local_u[:, 1] .= 1.0  # Left boundary
#     end
#     if rank == nranks - 1
#         local_u[:, Nlocal_j] .= 1.0  # Right boundary
#     end

#     global_diff = 0.0

#     # Main iteration loop
#     for iter in 1:MAX_ITER
#         # Perform one iteration of SOR
#         sor_iteration(local_u, Nlocal, rank, nranks)

#         # Exchange boundary rows with neighboring processes
#         exchange_boundaries(local_u, rank, nranks, Nlocal_i, Nlocal_j)

#         # Compute the local maximum residual
#         local_diff = 0.0
#         for i in 2:Nlocal_i-1
#             for j in 2:Nlocal_j-1
#                 local_diff = max(local_diff, abs(local_u[i, j] -
#                     0.25 * (local_u[i-1, j] + local_u[i+1, j] +
#                             local_u[i, j-1] + local_u[i, j+1])))
#             end
#         end

#         # Compute global maximum residual
#         global_diff = MPI.Allreduce(local_diff, MPI.MAX, MPI.COMM_WORLD)

#         # Output the residual at each iteration
#         if rank == 0
#             println("Iteration $iter, Residual: $global_diff")
#         end

#         # Check for convergence
#         if global_diff < TOL
#             if rank == 0
#                 println("Converged after $iter iterations")
#             end
#             break
#         end
#     end

#     MPI.Finalize()
# end

# # Run the main function
# main()


using MPI

# Function to perform SOR
function sor(A, b, x, omega, max_iter, tol)
    n = size(A, 1)
    for iter in 1:max_iter
        x_old = copy(x)
        for i in 1:n
            sum1 = A[i, 1:i-1] * x_old[1:i-1]
            sum2 = A[i, i+1:end] * x_old[i+1:end]
            x[i] = (1 - omega) * x_old[i] + (omega / A[i, i]) * (b[i] - sum1 - sum2)
        end
        # Check for convergence
        if norm(x - x_old) < tol
            break
        end
    end
    return x
end

function main()
    MPI.Init()

    # Get the rank and size
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    size = MPI.Comm_size(MPI.COMM_WORLD)

    # Matrix A and vector b setup (this should ideally be read from a file or generated)
    n = 4  # Number of equations (must be divisible by size)
    A = [4 1 2 0; 1 4 0 2; 2 0 4 1; 0 2 1 4]
    b = [1, 2, 3, 4]

    # Divide the rows of A and b among the processes
    rows_per_process = div(n, size)
    local_A = A[(rank * rows_per_process + 1):(rank + 1) * rows_per_process, :]
    local_b = b[(rank * rows_per_process + 1):(rank + 1) * rows_per_process]

    # Initial guess and parameters
    x = zeros(n)
    omega = 1.25
    max_iter = 1000
    tol = 1e-6

    # Perform SOR
    x = sor(local_A, local_b, x[(rank * rows_per_process + 1):(rank + 1) * rows_per_process], omega, max_iter, tol)

    # Gather results
    global_x = zeros(n)
    MPI.Allgather(x, global_x)

    if rank == 0
        println("Solution x: ", global_x)
    end

    MPI.Finalize()
end

main()
