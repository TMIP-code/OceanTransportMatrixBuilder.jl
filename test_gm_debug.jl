# Debug script to test GM implementation step by step
println("Starting GM implementation debug...")

try
    # Load the module explicitly
    using OceanTransportMatrixBuilder
    println("✓ Module loaded")
    
    # Check if functions are available
    println("Checking function availability:")
    println("  bolus_GM_velocity defined: ", @isdefined(bolus_GM_velocity))
    println("  buildTGM defined: ", @isdefined(buildTGM))
    
    # If not defined, try with module qualification
    if !@isdefined(bolus_GM_velocity)
        println("Trying with module qualification...")
        const bolus_GM_velocity = OceanTransportMatrixBuilder.bolus_GM_velocity
        const buildTGM = OceanTransportMatrixBuilder.buildTGM
        println("  bolus_GM_velocity defined: ", @isdefined(bolus_GM_velocity))
        println("  buildTGM defined: ", @isdefined(buildTGM))
    end
    
    # Create minimal test data
    nx, ny, nz = 3, 3, 2
    ρ = fill(1025.0, nx, ny, nz)
    ρ[2, 2, 1] = 1026.0
    
    # Create gridmetrics
    lon = range(0.0, 2.0, length=nx+1)
    lat = range(0.0, 2.0, length=ny+1)
    Z = [0.0, -50.0]
    gridtopology = (periodic_x=false, periodic_y=false)
    v3D = fill(1.0, nx, ny, nz)
    thkcello = fill(50.0, nx, ny, nz)
    edge_length_2D = (east=fill(1.0, nx, ny), north=fill(1.0, nx, ny))
    area2D = fill(1.0, nx, ny)
    distance_to_neighbour_2D = (west=fill(1.0, nx, ny), east=fill(1.0, nx, ny), 
                               south=fill(1.0, nx, ny), north=fill(1.0, nx, ny))
    
    gridmetrics = (; lon, lat, Z, gridtopology, v3D, thkcello, edge_length_2D, area2D, distance_to_neighbour_2D)
    indices = makeindices(v3D)
    
    println("✓ Test data created")
    
    # Test function calls
    println("Testing bolus_GM_velocity...")
    u_GM, v_GM = bolus_GM_velocity(ρ, gridmetrics, indices, κGM=600, maxslope=0.01)
    println("✓ bolus_GM_velocity succeeded")
    
    println("Testing buildTGM...")
    TGM = buildTGM(ρ, gridmetrics, indices, κGM=600, maxslope=0.01)
    println("✓ buildTGM succeeded")
    
    println("\n🎉 All GM functions working!")
    
catch e
    println("❌ Error: ", e)
    println("\nStack trace:")
    Base.show_backtrace(stdout, catch_backtrace())
end