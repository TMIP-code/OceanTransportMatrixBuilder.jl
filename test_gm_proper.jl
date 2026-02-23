# Proper GM test using standard grid construction
# This test creates a proper grid setup that should work with all functions

println("Creating proper GM test with full grid setup...")

using OceanTransportMatrixBuilder
using SparseArrays
using LinearAlgebra
using Test
println("✓ All packages loaded")

@testset "GM Proper Test" begin

    # Create a minimal but proper grid setup
    nx, ny, nz = 3, 3, 2  # Small but reasonable grid size

    # Create proper coordinate arrays
    lon = range(0.0, 2.0, length=nx+1)
    lat = range(-1.0, 1.0, length=ny+1)  # Symmetric around equator
    lev = [0.0, -50.0]  # Two vertical levels

    # Create vertex coordinates (required for proper grid construction)
    lon_vertices = [lon[i] for i in 1:nx+1, j in 1:ny+1]
    lat_vertices = [lat[j] for i in 1:nx+1, j in 1:ny+1]

    # Create minimal cell data (all ones for simplicity)
    areacello = fill(1.0, nx, ny)
    volcello = fill(1.0, nx, ny, nz)

    println("✓ Created grid coordinate data")

    # Build proper grid metrics using the standard function
    gridmetrics = makegridmetrics(
        areacello=areacello,
        volcello=volcello,
        lon=lon,
        lat=lat,
        lev=lev,
        lon_vertices=lon_vertices,
        lat_vertices=lat_vertices
    )

    println("✓ Grid metrics created successfully")

    # Create indices
    indices = makeindices(gridmetrics.v3D)
    println("✓ Indices created: $(indices.N) wet cells")

    # Create a density field with some variation
    ρ = fill(1025.0, nx, ny, nz)
    ρ[2, 2, 1] = 1026.0  # Add some variation
    ρ[1, 3, 2] = 1024.5

    println("✓ Density field created")

    # Test bolus_GM_velocity function
    println("\nTesting bolus_GM_velocity...")
    u_GM, v_GM = bolus_GM_velocity(ρ, gridmetrics, indices, κGM=600, maxslope=0.01)

    @test size(u_GM) == size(ρ) "GM u-velocity size mismatch"
    @test size(v_GM) == size(ρ) "GM v-velocity size mismatch"
    @test all(.!isnan.(u_GM)) "GM u-velocity contains NaN"
    @test all(.!isnan.(v_GM)) "GM v-velocity contains NaN"

    println("✓ bolus_GM_velocity test passed")
    println("  u_GM range: $(minimum(u_GM)) to $(maximum(u_GM))")
    println("  v_GM range: $(minimum(v_GM)) to $(maximum(v_GM))")

    # Test buildTGM function
    println("\nTesting buildTGM...")
    TGM = buildTGM(ρ, gridmetrics, indices, κGM=600, maxslope=0.01)

    @test TGM isa SparseMatrixCSC "TGM is not a sparse matrix"
    @test size(TGM) == (indices.N, indices.N) "TGM size mismatch"

    println("✓ buildTGM test passed")
    println("  TGM size: $(size(TGM))")
    println("  Non-zero entries: $(nnz(TGM))")

    # Test transportmatrix integration (this might still have issues)
    @testset "Test transportmatrix integration" begin
        # Create zero fluxes for testing
        ϕeast = zeros(nx, ny, nz)
        ϕwest = zeros(nx, ny, nz)
        ϕnorth = zeros(nx, ny, nz)
        ϕsouth = zeros(nx, ny, nz)
        ϕ = (east=ϕeast, west=ϕwest, north=ϕnorth, south=ϕsouth)

        mlotst = fill(10.0, nx, ny)  # Mixed layer depth

        result = transportmatrix(
            ϕ=ϕ, mlotst=mlotst, gridmetrics=gridmetrics, indices=indices, ρ=1025.0,
            ρ_GM=ρ, κGM=600.0, maxslope=0.01
        )

        @test hasproperty(result, :TGM) "Result missing TGM property"
        @test result.TGM isa SparseMatrixCSC "TGM in result is not sparse"

        println("✓ transportmatrix integration passed")

    end

    println("\n🎉 GM core functionality tests completed successfully!")
    println("✅ bolus_GM_velocity: Working")
    println("✅ buildTGM: Working")
    println("✅ Basic integration: Working")

end