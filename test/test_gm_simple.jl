# Simple test for GM density-only implementation
# This test verifies that the basic functionality works without requiring full grid setup

using Test
using OceanTransportMatrixBuilder
using SparseArrays
using LinearAlgebra

@testset "GM density-only basic functionality" begin
    # Create minimal test data
    nx, ny, nz = 3, 3, 2
    
    # Create a simple density field with some variation
    # Use interior cells only to avoid boundary issues
    ρ = fill(1025.0, nx, ny, nz)  # Base density 1025 kg/m³
    ρ[2, 2, 1] = 1026.0  # Add variation in interior cell
    ρ[2, 2, 2] = 1024.5  # Add variation in same column
    
    # Create minimal gridmetrics (simplified for testing)
    # In a real scenario, these would be properly constructed
    lon = range(0.0, 2.0, length=nx+1)
    lat = range(0.0, 2.0, length=ny+1)
    Z = [0.0, -50.0]  # Simple 2-level vertical grid
    
    # Create a simple grid topology (this would normally be computed)
    # For testing, we'll create a minimal grid topology that provides the needed functions
    struct SimpleGridTopology <: OceanTransportMatrixBuilder.AbstractGridTopology
        nx::Int
        ny::Int
        nz::Int
    end
    
    # Define the required grid navigation functions for our simple grid
    OceanTransportMatrixBuilder.i₊₁(C, g::SimpleGridTopology) = C.I[1] < g.nx ? C + CartesianIndex(1, 0, 0) : nothing
    OceanTransportMatrixBuilder.i₋₁(C, g::SimpleGridTopology) = C.I[1] > 1 ? C + CartesianIndex(-1, 0, 0) : nothing
    OceanTransportMatrixBuilder.j₊₁(C, g::SimpleGridTopology) = C.I[2] < g.ny ? C + CartesianIndex(0, 1, 0) : nothing
    OceanTransportMatrixBuilder.j₋₁(C, g::SimpleGridTopology) = C.I[2] > 1 ? C + CartesianIndex(0, -1, 0) : nothing
    OceanTransportMatrixBuilder.k₊₁(C, g::SimpleGridTopology) = C.I[3] < g.nz ? C + CartesianIndex(0, 0, 1) : nothing
    OceanTransportMatrixBuilder.k₋₁(C, g::SimpleGridTopology) = C.I[3] > 1 ? C + CartesianIndex(0, 0, -1) : nothing
    
    gridtopology = SimpleGridTopology(nx, ny, nz)
    
    # Create volume and thickness fields
    v3D = fill(1.0, nx, ny, nz)  # Volume
    thkcello = fill(50.0, nx, ny, nz)  # Layer thickness
    edge_length_2D = (east=fill(1.0, nx, ny), north=fill(1.0, nx, ny))
    area2D = fill(1.0, nx, ny)
    distance_to_neighbour_2D = (west=fill(1.0, nx, ny), east=fill(1.0, nx, ny), 
                               south=fill(1.0, nx, ny), north=fill(1.0, nx, ny))
    
    gridmetrics = (; lon, lat, Z, gridtopology, v3D, thkcello, edge_length_2D, area2D, distance_to_neighbour_2D)
    
    # Create indices
    indices = makeindices(v3D)
    
    @testset "bolus_GM_velocity function" begin
        # Test that bolus_GM_velocity runs without error
        u_GM, v_GM = bolus_GM_velocity(ρ, gridmetrics, indices, κGM=600, maxslope=0.01)
        
        # GM u-velocity should have same size as density field
        @test size(u_GM) == size(ρ)
        # GM v-velocity should have same size as density field
        @test size(v_GM) == size(ρ)
        # GM u-velocity should not contain NaN values
        @test all(.!isnan.(u_GM))
        # GM v-velocity should not contain NaN values
        @test all(.!isnan.(v_GM))
        
        # Test that velocities are not all zero (should respond to density variation)
        # GM should produce non-zero velocities
        @test any(!isapprox(u_GM[i], 0.0, atol=1e-10) for i in eachindex(u_GM))
        @test any(!isapprox(v_GM[i], 0.0, atol=1e-10) for i in eachindex(v_GM))
    end
    
    @testset "buildTGM function" begin
        # Test that buildTGM runs without error
        TGM = buildTGM(ρ, gridmetrics, indices, κGM=600, maxslope=0.01)
        
        # TGM should be a sparse matrix
        @test TGM isa SparseMatrixCSC
        # TGM should be square matrix of correct size
        @test size(TGM) == (indices.N, indices.N)
        
        # Test that TGM is not a zero matrix when density varies
        # TGM should have non-zero entries when density varies
        @test nnz(TGM) > 0
        
        # Test that TGM has reasonable properties
        # TGM should not contain NaN values
        @test all(.!isnan.(TGM.nzval))
        # TGM should not contain Inf values
        @test all(.!isinf.(TGM.nzval))
    end
    
    @testset "transportmatrix with GM" begin
        # Create minimal flux data for testing
        # In a real scenario, these would be properly computed from velocity data
        ϕeast = zeros(nx, ny, nz)
        ϕwest = zeros(nx, ny, nz)
        ϕnorth = zeros(nx, ny, nz)
        ϕsouth = zeros(nx, ny, nz)
        ϕ = (east=ϕeast, west=ϕwest, north=ϕnorth, south=ϕsouth)
        
        # Test transportmatrix with GM
        mlotst = fill(10.0, nx, ny)  # Mixed layer depth
        
        # This should work without error
        result = transportmatrix(
            ϕ=ϕ, mlotst=mlotst, gridmetrics=gridmetrics, indices=indices, ρ=1025.0,
            ρ_GM=ρ, κGM=600.0, maxslope=0.01
        )
        
        # Result should include TGM operator
        @test hasproperty(result, :TGM)
        # TGM in result should be sparse matrix
        @test result.TGM isa SparseMatrixCSC
        # TGM should have correct size
        @test size(result.TGM) == (indices.N, indices.N)
        
        # Test that total transport matrix includes GM contribution
        T_total = result.T
        @test T_total ≈ result.Tadv + result.TκH + result.TκVML + result.TκVdeep + result.TGM
    end
end

@info "GM density-only tests completed successfully!"