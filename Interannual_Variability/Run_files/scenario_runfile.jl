##################### Sample Runfile for GenX runs #####################
using GenX
using JuMP
using DataFrames
using CSV

file_path = @__DIR__
hm_dir = dirname(file_path)

##### Folder structure - one major folder with three sub-folders: Scenarios, Runfiles, and Results

## Get the scenario length and set directory name
scenario_length = 10
scen_directory = scenario_length * "_length_scenarios"
res_dir = scenario_length * "_results"

## Get the path to the inputs folder and make the corresponding results folder
run_path = joinpath(hm_dir, "Scenarios", scen_directory)
res_path_gf = joinpath(hm_dir, "Results", res_dir * "_gf")
res_path_bf = joinpath(hm_dir, "Results", res_dir * "_bf")

mkdir(res_path_gf)
mkdir(res_path_bf) 

## Grab the run_helpers.jl from the directory
include(joinpath(file_path, "run_helpers.jl"))

## Initialize the settings (same for all of one type of run)
settings_path = get_settings_path(run_path)
genx_settings = get_settings_path(run_path, "genx_settings.yml")
mysetup = configure_settings(genx_settings)

## Get a vector of strings that has the names of the folders within the inputs folder
scen_folders = readdir(run_path)

## The list of emissions constraints (in gCO2/kWh)
emiss_lim_list = [12.0]

## Make a folder for each emissions limit
for emiss_cap in emiss_lim_list
    ## Make the folder
    emiss_folder_gf = joinpath(res_path_gf, string(emiss_cap) * "_limit")
    emiss_folder_bf = joinpath(res_path_bf, string(emiss_cap) * "_limit")
    mkdir(emiss_folder_gf)
    mkdir(emiss_folder_bf)
end 

## Loop through each folder and run the GenX case
for folder in scen_folders
    ## Set inputs path
    inputs_path = joinpath(run_path, folder)

    ## Cluster time series inputs if necessary and if specified by the user
    TDRpath = joinpath(inputs_path, mysetup["TimeDomainReductionFolder"])
    if mysetup["TimeDomainReduction"] == 1
        if !time_domain_reduced_files_exist(TDRpath)
            println("Clustering Time Series Data (Grouped)...")
            cluster_inputs(inputs_path, settings_path, mysetup)
        else
            println("Time Series Data Already Clustered.")
        end
    end

    ## Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)

    # Turn this setting on if you run into numerical stability issues
    # set_optimizer_attribute(OPTIMIZER, "BarHomogeneous", 1)

    #### Running a case

    ## Load inputs
    println("Loading Inputs")
    myinputs = load_inputs(mysetup, inputs_path)

    ## Set the load-based emissions cap
    mysetup["CO2Cap"] = 2

    ## Set the number of CO2CapPeriods to the length of the scenario
    mysetup["CO2CapPeriods"] = scenario_length

    ## Set the scale factor
    scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    ## Run the case for each emissions limit you have
    for emiss_limit in emiss_lim_list
        ## Set the output dir_name
        output_dirname = string(folder) * "_" * string(emiss_limit) * "_limit"
        output_path = joinpath(res_path_gf, string(emiss_limit) * "_limit", output_dirname) 

        mkdir(output_path)

        # Hard-coded to put all emissions in New Hampshire, but CO2 Cap is set to be system-wide
        myinputs["dfMaxCO2Rate"][2] = emiss_lim / scale_factor ./ 1e3

        if isfile(joinpath(outputs_path, "costs.csv"))
            println("Skipping Case for emiss limit = " * string(emiss_lim) * " because it already exists.")
            continue
        end

        ## Generate model
        println("Generating the Optimization Model")
        EP = generate_model(mysetup, myinputs, OPTIMIZER)

        ########################
        #### Add any additional constraints
        HYDRO_RES = myinputs["HYDRO_RES"]
        dfGen = myinputs["dfGen"]

        # Empty arrays for indexing
        jan1_idxs = Int[]
        may1_idxs = Int[]

        # year indexing
        for year_num in 1:scenario_length
            # Calculate the index for the beginning of the years
            start_year = (year_num-1) * 8760 + 1
            push!(jan1_idxs, start_year)

            # Calculate the index for the middle of the years
            mid_year = (year_num-1) * 8760 + 2879
            push!(may1_idxs, mid_year)
        end

        ## Hydro storage == 0.70 * Existing Capacity at the start of the year
        @constraint(EP, cHydroJan[y in HYDRO_RES, jan1_idx in jan1_idxs], EP[:vS_HYDRO][y, jan1_idx]  .== 0.70 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio])
        
        ## Hydro storage <= 0.55 * Existing Capacity at start of May 1st 
        @constraint(EP, cHydroSpring[y in HYDRO_RES, may1_idx in may1_idxs], EP[:vS_HYDRO][y, may1_idx] .<= 0.55 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio])
        
        ## Maine -> Quebec transmission limited to 2170MWe.
        # The line is defined as Quebec -> Maine in Network.csv, so these flows will be negative
        # Make sure to correc the line index if the order is changed in Network.csv
        @constraint(EP, cMaine2Quebec[t=1:myinputs["T"]], EP[:vFLOW][2, t] >= -170.0)

        ## Solve model
        println("Solving Model")
        EP, solve_time = solve_model(EP, mysetup)
        myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

        ## Run MGA if the MGA flag is set to 1 else only save the least cost solution
        println("Writing Output")

        ## Write only the outputs of interest: costs, capacity, and capacity_factor
        write_capacity(output_path, myinputs, mysetup, EP)
        println("Capacity written for greenfield")
        write_capacityfactor(path, myinputs, mysetup, EP)
        println("Capacity factor written for greenfield")
        write_costs(path, myinputs, mysetup, EP)
        println("Costs written for greenfield")

        ################################################################################################ BROWNFIELD CASES

        #### Run the corresponding brownfield case ####
        ## Grab the run_helpers.jl from the directory
        include(joinpath(file_path, "run_helpers.jl"))

        bf_inputs_path = joinpath(home_dir, "brownfield_data")    ## This is where the inputs for the brownfield are stored

        ## Initialize the settings (same for all of one type of run)
        bf_settings_path = get_settings_path(bf_inputs_path)
        bf_genx_settings = get_settings_path(bf_inputs_path, "genx_settings.yml")
        mysetup_bf = configure_settings(bf_genx_settings)
        
        ## Cluster time series inputs if necessary and if specified by the user
        TDRpath = joinpath(bf_inputs_path, mysetup_bf["TimeDomainReductionFolder"])
        if mysetup["TimeDomainReduction"] == 1
            if !time_domain_reduced_files_exist(TDRpath)
                println("Clustering Time Series Data for Brownfield (Grouped)...")
                cluster_inputs(bf_inputs_path, bf_settings_path, mysetup_bf)
            else
                println("Time Series Data Already Clustered.")
            end
        end

        ## Configure solver
        println("Configuring Solver for Brownfield")
        OPTIMIZER_bf = configure_solver(mysetup_bf["Solver"], bf_settings_path)

        # Turn this setting on if you run into numerical stability issues
        # set_optimizer_attribute(OPTIMIZER, "BarHomogeneous", 1)

        #### Running a case

        ## Load inputs
        println("Loading Inputs")
        myinputs_bf = load_inputs(mysetup_bf, bf_inputs_path)

        ## Grab the corresponding output capacity and use that to run a brownfield
        cap_out = CSV.read(joinpath(output_path, "capacity.csv"), DataFrame)

        ## Grab the power and energy capacities
        yearly_cap = cap_out[1:end-1, "EndCap"]
        yearly_mwh = cap_out[1:end-1, "EndEnergyCap"]

        ## Set them in the inputs
        myinputs_bf["dfGen"][!, "Existing_Cap_MW"] = yearly_cap
        myinputs_bf["dfGen"][!, "Existing_Cap_MWh"] = yearly_mwh

        ## Set the load-based emissions cap
        mysetup_bf["CO2Cap"] = 2   

        ## Set the number of CO2CapPeriods to the length of the scenario
        mysetup_bf["CO2CapPeriods"] = 23  

        ## Set the scale factor
        scale_factor_bf = mysetup_bf["ParameterScale"] == 1 ? ModelScalingFactor : 1

        # Hard-coded to put all emissions in New Hampshire, but CO2 Cap is set to be system-wide
        myinputs_bf["dfMaxCO2Rate"][2] = emiss_lim / scale_factor_bf ./ 1e3

        # Set the output path for bf 
        output_path_bf = joinpath(res_path_bf, string(emiss_limit) * "_limit", output_dirname)
        mkdir(output_path_bf)

       if isfile(joinpath(output_path_bf, "costs.csv"))
           println("Skipping Case for emiss limit = " * string(emiss_lim) * " because it already exists.")
           continue
       end

       ## Generate model
       println("Generating the Optimization Model")
       EP_bf = generate_model(mysetup_bf, myinputs_bf, OPTIMIZER_bf)

       ########################
       #### Add any additional constraints
       HYDRO_RES_bf = myinputs_bf["HYDRO_RES"]
       dfGen_bf = myinputs_bf["dfGen"]

       # Empty arrays for indexing
       jan1_idxs_bf = Int[]
       may1_idxs_bf = Int[]

       # year indexing
       for year_num in 1:23
           # Calculate the index for the beginning of the years
           start_year = (year_num-1) * 8760 + 1
           push!(jan1_idxs_bf, start_year)

           # Calculate the index for the middle of the years
           mid_year = (year_num-1) * 8760 + 2879
           push!(may1_idxs_bf, mid_year)
       end

       ## Hydro storage == 0.70 * Existing Capacity at the start of the year
       @constraint(EP_bf, cHydroJan[y in HYDRO_RES_bf, jan1_idx_bf in jan1_idxs_bf], EP_bf[:vS_HYDRO][y, jan1_idx_bf]  .== 0.70 .* EP_bf[:eTotalCap][y] .* dfGen_bf[y,:Hydro_Energy_to_Power_Ratio])
       
       ## Hydro storage <= 0.55 * Existing Capacity at start of May 1st 
       @constraint(EP_bf, cHydroSpring[y in HYDRO_RES_bf, may1_idx_bf in may1_idxs_bf], EP_bf[:vS_HYDRO][y, may1_idx_bf] .<= 0.55 .* EP_bf[:eTotalCap][y] .* dfGen_bf[y,:Hydro_Energy_to_Power_Ratio])
       
       ## Maine -> Quebec transmission limited to 2170MWe.
       # The line is defined as Quebec -> Maine in Network.csv, so these flows will be negative
       # Make sure to correc the line index if the order is changed in Network.csv
       @constraint(EP_bf, cMaine2Quebec[t=1:myinputs["T"]], EP_bf[:vFLOW][2, t] >= -170.0)

       ## Solve model
       println("Solving Model")
       EP_bf, solve_time = solve_model(EP_bf, mysetup_bf)
       myinputs_bf["solve_time"] = solve_time # Store the model solve time in myinputs

       ## Run MGA if the MGA flag is set to 1 else only save the least cost solution
       println("Writing Output")

       ## Write only the outputs of interest: costs, capacity, and capacity_factor
       write_capacity(output_path_bf, myinputs_bf, mysetup_bf, EP_bf)
       println("Capacity written for brownfield")
       write_capacityfactor(output_path_bf, myinputs_bf, mysetup_bf, EP_bf)
        println("Capacity factor written for brownfield")
       write_costs(output_path_bf, myinputs_bf, mysetup_bf, EP_bf)
       println("Costs written for brownfield")
    end
end


################################
## Add folder for brownfield data
## Add folder for brownfield results
## Add Nuclear to Generators_data.csv
## Finish making Generators_variability.csv
        