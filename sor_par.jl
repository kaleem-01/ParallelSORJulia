using MPI

function sor_par!(A, stopdiff, maxiters, comm)
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    ni, nj = size(A)

    # Ensure `ni` is a multiple of `nranks`
    # if mod(ni, nranks) != 0
    #     println("Error: `ni` must be a multiple of `nranks`.")
    #     MPI.Abort(comm, 1)
    # end

    n_own = div(ni, nranks)
    
    iter_count = 0
    ω = 0.25
    diff = zero(eltype(A))
    for iter in 1:maxiters
        reqs = MPI.Request[]
        if nranks == 1
            for i in 2:(ni-1)
                for j in 2:(nj-1)
                    Aij_old = A[i,j]
                    Aij_new = ω*(A[i-1,j]+A[i+1,j]+A[i,j-1]+A[i,j+1])
                    diff_Aij = abs(Aij_new - Aij_old)
                    diff = max(diff, diff_Aij)
                    A[i, j] = Aij_new
                end
            end
            iter_count += 1
            if diff < stopdiff
                break
            end
        # Communicate boundary rows with neighbors
        elseif rank == 0
            # Send last row to next rank and receive from it
            req = MPI.Isend(view(A, n_own, 2:nj-1), comm, dest=1, tag=0)
            push!(reqs, req)
            req = MPI.Irecv!(view(A, n_own + 1, 2:nj-1), comm, source=1, tag=0)
            push!(reqs, req)
        elseif rank == nranks - 1
            # Send first row to previous rank and receive from it
            req = MPI.Isend(view(A, rank * n_own + 1, 2:nj-1), comm, dest=rank - 1, tag=0)
            push!(reqs, req)
            req = MPI.Irecv!(view(A, rank * n_own, 2:nj-1), comm, source=rank - 1, tag=0)
            push!(reqs, req)
        else
            # Middle ranks: communicate with both top and bottom neighbors
            top_neigh = rank - 1
            bottom_neigh = rank + 1
            
            req = MPI.Isend(view(A, rank * n_own + 1, 2:nj-1), comm, dest=top_neigh, tag=0)
            push!(reqs, req)
            req = MPI.Isend(view(A, (rank + 1) * n_own, 2:nj-1), comm, dest=bottom_neigh, tag=0)
            push!(reqs, req)
            
            req = MPI.Irecv!(view(A, rank * n_own, 2:nj-1), comm, source=top_neigh, tag=0)
            push!(reqs, req)
            req = MPI.Irecv!(view(A, (rank + 1) * n_own + 1, 2:nj-1), comm, source=bottom_neigh, tag=0)
            push!(reqs, req)
        end

        # Wait for all communications to complete
        MPI.Waitall(reqs)
        
        # Perform computation on inner rows
        local_diff = zero(eltype(A))
        for i in (rank * n_own + 2):(rank + 1) * n_own - 1
            for j in 2:(nj - 1)
                Aij_old = A[i, j]
                Aij_new = ω * (A[i-1, j] + A[i+1, j] + A[i, j-1] + A[i, j+1])
                diff_Aij = abs(Aij_new - Aij_old)
                local_diff = max(local_diff, diff_Aij)
                A[i, j] = Aij_new
            end
        end

        # Reduce the maximum `diff` value across all ranks
        global_diff = MPI.Allreduce(local_diff, MPI.MAX, comm)
        
        # Update iteration count and check for convergence
        iter_count += 1
        if global_diff < stopdiff
            break
        end
    end
    return stopdiff, iter_count     
end

# Initialize MPI and test the function
# MPI.Init()
# comm = MPI.COMM_WORLD
# rank = MPI.Comm_rank(comm)

# A = rand(100, 100)
# stopdiff = 1e-6
# maxiters = 1000

# # Call the parallel SOR function
# result_stopdiff, result_iters = sor_par!(A, stopdiff, maxiters, comm)

# if rank == 0
#     println("Converged after $result_iters iterations with final diff: $result_stopdiff")
# end

# MPI.Finalize()

