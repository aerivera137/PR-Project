# Exterior functions for hurricane
# -> Get wind speed
function get_windspeed()
    # Select wind speed of hurricane
    d_wind = TriangularDist(75,158,150) # lower limit, upper limit, most are Category 2
    ws = rand(d_wind)
    return ws
end
# -> Generate hurricane category + period (hours)
function generate_hurricane(time_periods)    
    # Duration of hurricane season (when a hurricane can occur)
    JUNE_1_DAY = 152 # hour: 2921
    NOV_30_DAY = 334 # hour: 8030
    # Select a day to have a hurricane (higher probability of occurring in September)
    hours_per_day = 24
    FIRST_DAY_SEPT = 244
    d_day = TriangularDist(JUNE_1_DAY,NOV_30_DAY,FIRST_DAY_SEPT) # upper limit, lower limit, most in September
    HURRICANE_DAY_START = rand(d_day) 
    HURRICANE_DAY_START = convert(Int64,round(HURRICANE_DAY_START))
    HURRICANE_HOUR_START = HURRICANE_DAY_START * hours_per_day
    # Set duration of outage 
    P_HURRICANE = HURRICANE_HOUR_START:(HURRICANE_HOUR_START+48)
    # Get wind speed
    wind_speed = get_windspeed()
    # Discretize wind speed into a category of hurricane (where 2/3 and 4/5 are combined)
    hurricane_category_boundaries = [75,96,111,130,157]
    category = searchsortedfirst(hurricane_category_boundaries, wind_speed) - 1
    return category, P_HURRICANE
end
# Processing inputs
# -> Microgrids
function process_MG_inputs(generators,zone)
    # Define sets of microgrids in each zone (1=East,2=North,3=South,4=West)
    MG_Zone = Vector{Int64}()

    for (i, row) in enumerate(eachrow(generators))
        if row[:Zone] == zone && (contains(row[:technology], "Microgrid_Solar_Moderate") || contains(row[:technology], "Microgrid_Diesel_Moderate"))
            push!(MG_Zone, i)
        end
    end

    return MG_Zone
end
"# -> Transmission lines
function process_L_inputs(all_inputs)
    # Define sets of lines

    return NVL_set, VL_set, BL_set
end"
# -> Failure rates for each generator in each category of hurricane
function get_failure_rate(generators,gen,cat)
    column = "Cat_" * string(cat) * "_Failure_Rate"
    return generators[gen, column]
end
# -> Limit max microgrid power output when it's protecting a substation
function limit_mg_output!(mg,generators,inputs,P_HURRICANE,zone,mg_count)
    # Get its capacity in MW
    mg_capacity = generators[mg,"Existing_Cap_MW"]
    # Get the demand profile during the hurricane in MW
    load_hurricane = Dict{Int,Any}()
    for hour in P_HURRICANE
        load_hurricane[hour] = inputs["pD"][hour,zone]
    end
    # Set the correct indices of inputs["pP_Max"] to the ratio of these, normalize by number of microgrids available in this zone
    for hour in P_HURRICANE
        if load_hurricane[hour] != 0
            inputs["pP_Max"][mg,hour] = (load_hurricane[hour] / mg_capacity) / mg_count
        end
    end
end
# -> Modify the load in a zone that doesn't have a microgrid
function modify_load!(inputs,original_load,zone,P_HURRICANE,num_substations)
    # Modify load at each hour 
    for hour in P_HURRICANE
        # Un-serve demand in this zone for these hours 
        inputs["pD"][hour,zone] -= (1/num_substations)*original_load[hour,zone]
    end
end
# Main function: running hurricane simulation
function hurricane_sim!(model::Model, inputs::Dict)
    ####################################################################################################################################
    ############################################################## SET-UP ##############################################################
    ####################################################################################################################################
    println("Starting hurricane sim")
    #= NEED TO DEFINE: 
    - T (hours)
    - Z (zones)
    - G (all generators), UC (thermal gens), VRE (VREs); MG_East, MG_North, MG_South, MG_West (microgrids in each zone)
    - L (all lines), NVL (non-vulnerable lines), VL (vulnerable lines), BL (buried lines)
    =#

    T = 1:inputs["T"]
    G = 1:inputs["G"]
    UC = inputs["COMMIT"] 
    VRE = inputs["VRE"] 
    dfGen = inputs["dfGen"]
    MG_East = process_MG_inputs(dfGen,1)
    MG_North = process_MG_inputs(dfGen,2)
    MG_South = process_MG_inputs(dfGen,3)
    MG_West = process_MG_inputs(dfGen,4)
    Z = 1:inputs["Z"]
    L = 1: inputs["L"]
    "NVL, VL, BL = process_L_inputs(inputs)"

    # Copy the original demand profile as a reference
    original_load = copy(inputs["pD"])

    # Define substation failure rates
    Substation_FailRates = DataFrame(Category=[1,2,3,4,5], Failure_Rate=[0.000000626,0.002837392,0.002837392,0.269302475,0.269302475])

    # Define default transmission line failure rates
    TranLine_FailRates = DataFrame(Category=[1,2,3,4,5], Failure_Rate=[0.001005387,0.036422691,0.036422691,0.27497107,0.27497107])

    # Define vulnerable transmission line failure rates (20% more vulnerable than default)
    VulnLine_FailRates = DataFrame(Category=[1,2,3,4,5], Failure_Rate=[0.001206464,0.043707229,0.043707229,0.329965284,0.329965284])

    # Define buried transmission line failure rates
    BurLine_FailRates = DataFrame(Category=[1,2,3,4,5], Failure_Rate=[0,0,0,0.15841923,0.158419238])

    # Define an array of outage values for each UC generator (0 for no outage, 1 for outage)
    outages = zeros(length(T), length(G))  

    "# Define an array of damage values for each transmission line (1 for no damage, 0 for damage)
    line_availability = ones(length(T), length(L))"

    # Check if substations are secure, if not, some demand is unmet in that area
    function check_substations(z, MG_Zone, generators, original_load, P_HURRICANE, category)
        println("Checking substations")
        # Define peak demand in the zone (MW)
        peak_demand_in_zone = maximum(inputs["pD"][:, z])*1000
        println("Peak demand in zone:")
        println(string(peak_demand_in_zone))
        # Define number of substations in the zone based on peak demand, where all substations (except remainder) have 20 MW
        num_substations = trunc(Int,peak_demand_in_zone/20) + convert(Int,rem(peak_demand_in_zone,20)>0)
        println("Number of substations:")
        println(string(num_substations))
        # Record capacity of remaining substation
        cap_remainder_substation = rem(peak_demand_in_zone,20)
        # Define number of microgrids in Eastern zone (maximum is 2 = 1 solar + 1 diesel)
        num_MG = length(MG_Zone)
        # Instantiate the number + capacity of available/working microgrids as 0
        num_MG_Available = 0
        cap_MG_Available = 0
        # Create array of all generators where a value of '1' means the microgrid is available, '0' means unavailable
        array_MG_Available = zeros(length(G))

        # Check if there are microgrids in the specified zone 
        if (num_MG != 0)
            # Check pairs of microgrids + substations (total of 10 substations)
            for mg in MG_Zone
                # Check if microgrid has failed
                if rand(Bernoulli(get_failure_rate(generators,mg,category)))==0
                    # Record number of working/available microgrids
                    println("There is a microgrid protecting a substation")
                    num_MG_Available += 1
                    cap_MG_Available += dfGen[mg,"Existing_Cap_MW"]
                    array_MG_Available[mg] = 1
                end
            end
        end

        # Record indices of microgrids that are available
        indices_MG_Available = findall(x -> x == 1, array_MG_Available)   

        # Limit microgrid capacity during the hurricane period if it's protecting a substation
        for mg in indices_MG_Available
            println("Limiting microgrid output")
            limit_mg_output!(mg,dfGen,inputs,P_HURRICANE,z,num_MG_Available)
        end

        # Check substations that don't have microgrids
        num_remaining_substations = num_substations - num_MG_Available
        println("Number of remaining substations:")
        println(string(num_remaining_substations))

        # Iterate through remaining substations
        for i in 1:num_remaining_substations
            # Check if remaining substation has failed
            if rand(Bernoulli(Substation_FailRates.Failure_Rate[category]))
                # If substation has failed, remove portion of demand for substation
                println("Original load:")
                println(inputs["pD"][P_HURRICANE, z])
                println("A substation has failed - modifying load")
                modify_load!(inputs, original_load, z, P_HURRICANE, num_substations)
                println(inputs["pD"][P_HURRICANE, z])
            end
        end

        #return sum(original_load[:, z])-sum(inputs["pD"][:, z])
        return original_load[:, z] - inputs["pD"][:, z]
    end

    # Function to apply damage to system components for a given category of hurricane
    function damage_sim(category,P_HURRICANE)
        println(string(first(P_HURRICANE)))
        println(string(last(P_HURRICANE)))
        println(string(category))
        # Change variability for variable renewable generators to 0 if there is an outage
        for g in VRE
            # Bernoulli draw returns 1 for generators with higher probability of outage
            if rand(Bernoulli(get_failure_rate(dfGen,g,category))) 
                # Force variability to zero for VREs that are out for the duration of the outage
                println("A VRE has failed - modifying variability")
                START = first(P_HURRICANE)
                END = last(P_HURRICANE)
                restoration_time = dfGen[g,"Hours_to_Restore"]
                outage_time = START:(END+restoration_time)
                inputs["pP_Max"][g,outage_time] .= 0
            end
        end

        # Update outage value for thermal UC generators to 1 if there is an outage
        for g in UC
            # Bernoulli draw returns 1 for generators with higher probability of outage
            if rand(Bernoulli(get_failure_rate(dfGen,g,category)))
                # Force generator to shut off for the duration of the outage
                println("A generator has failed - modifying commitment state")
                START = first(P_HURRICANE)
                END = last(P_HURRICANE)
                restoration_time = dfGen[g,"Hours_to_Restore"]
                outage_time = START:(END+restoration_time)
                outages[outage_time,g] .= 1
            end
        end
        
        # Check all substations for failure + protect demand if microgrids are available
        lost_load_z1 = check_substations(1,MG_East,dfGen,original_load,P_HURRICANE,category)
        lost_load_z2 = check_substations(2,MG_North,dfGen,original_load,P_HURRICANE,category)
        lost_load_z3 = check_substations(3,MG_South,dfGen,original_load,P_HURRICANE,category)
        lost_load_z4 = check_substations(4,MG_West,dfGen,original_load,P_HURRICANE,category)

        dfNSE = DataFrame(Lost_Load_Zone_1 = lost_load_z1, 
                          Lost_Load_Zone_2 = lost_load_z2, 
                          Lost_Load_Zone_3 = lost_load_z3,
                          Lost_Load_Zone_4 = lost_load_z4)

        "#Check for transmission line failure
        # Default/non-vulnerable lines
        for l in NVL
            # Check for failure
            if rand(Bernoulli(TranLine_FailRates.Failure_Rate[category])) 
                # If failed, cut off power flow in this line during the hurricane
                line_availability[P_HURRICANE,l] .= 0
            end
        end
        # Vulnerable lines
        for l in VL  
            # Check for failure
            if rand(Bernoulli(VulnLine_FailRates.Failure_Rate[category])) 
                # If failed, cut off power flow in this line during the hurricane
                line_availability[P_HURRICANE,l] .= 0
            end
        end
        # Buried lines
        for l in BL
            # Check for failure (these lines are in the generators data set)
            if rand(Bernoulli(get_failure_rate(generators,l,category)))
                # If failed, cut off power flow in this line during the hurricane
                line_availability[P_HURRICANE,l] .= 0
            end
        end"

        # Define an array with just the indices of UC generators that have outages
        #outage_indices = findall(x -> x == 1, outages)   
        return outages, dfNSE
    end
    ####################################################################################################################################
    ############################################################ SIMULATION ############################################################
    ####################################################################################################################################
    # Generate hurricane period and category
    category, P_HURRICANE = generate_hurricane(T)

    # Call damage function (main function)
    outages, dfNSE = damage_sim(category, P_HURRICANE)

    println("Returning outputs")
    return outages, dfNSE
end