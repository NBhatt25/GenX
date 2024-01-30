##########
# This Julia script reads in the scenario_list.csv file and uses it to create
# new folders that have the necessary files to run a GenX case for that scenario
##########

using CSV, DataFrames

#### Put this file in the main directory
main_dir = dirname(@__DIR__)

## Set the directory where the base files are located
base_files = joinpath(main_dir, "base_scenario_data")

## Set the directory where input data is located
input_data = joinpath(main_dir, "summary_data")

## Scenario files
scenario_dir = joinpath(main_dir, "Scenarios")

## First we want to open up the scenario list
scenario_list = CSV.read(joinpath(input_data, "scenario_list.csv"), DataFrame)

## Next, we want to open up the files that have the data we are interested in (load, solar, and wind)
load_data_csv = CSV.read(joinpath(input_data, "yearly_historical_load.csv"), DataFrame)
onsw_data_csv = CSV.read(joinpath(input_data, "yearly_historical_onsw.csv"), DataFrame)
ofsw_data_csv = CSV.read(joinpath(input_data, "yearly_historical_ofsw.csv"), DataFrame)
solar_data_csv = CSV.read(joinpath(input_data, "yearly_historical_solar.csv"), DataFrame)

## Then, we want to loop through each row of the scenario list and convert each row to a list
set_of_scenarios = []

for row in eachrow(scenario_list)
    push!(set_of_scenarios, collect(row))
end

## Create a new folder based on the length of the scenario that you have
scenario_length = length(set_of_scenarios[1])
dir_name = string(scenario_length) * "_length_scenarios"
dir_location = joinpath(scenario_dir, dir_name)

## Make the folder in the chosen directory
mkdir(dir_location)

### Add the Settings folder (ONLY NEED ONE FOR ALL SCENARIOS OF ONE LENGTH)
mkdir(joinpath(dir_location, "Settings"))

## Path for settings files
base_settings = joinpath(base_files, "Settings")
dest_settings = joinpath(dir_location, "Settings")

## Add the yml files to the Settings folder
cp((joinpath(base_settings, "genx_settings.yml")), (joinpath(dest_settings, "genx_settings.yml")))
cp((joinpath(base_settings, "gurobi_settings.yml")), (joinpath(dest_settings, "gurobi_settings.yml")))
cp((joinpath(base_settings, "time_domain_reduction_settings.yml")), (joinpath(dest_settings, "time_domain_reduction_settings.yml")))


## Now, we want to loop through each scenario and stitch together the data for those files
for scenario in set_of_scenarios
    ## Create a new folder for each scenario
    if scenario_length == 23
        scen_dirname = "scenario_" * "23yr"
    else
        scen_dirname = "scenario_" * join(scenario, "_")    
    end

    ## Make this directory in the folder that you created before looping
    scen_dirname_path = joinpath(dir_location, scen_dirname)
    mkdir(scen_dirname_path)
    
    ## Create an empty list to store the yearly data for each file
    load_data = []
    onsw_data = []
    ofsw_data = []
    solar_data = []

    for item in scenario
        ## Grab the corresponding column for each item in the scenario
        item_load = getproperty(load_data_csv, Symbol(item))
        item_onsw = getproperty(onsw_data_csv, Symbol(item + 2000))
        item_ofsw = getproperty(ofsw_data_csv, Symbol(item + 2000))
        item_solar = getproperty(solar_data_csv, Symbol(item + 2000))

        ## Append the data to the corresponding list
        push!(load_data, item_load)
        push!(onsw_data, item_onsw)
        push!(ofsw_data, item_ofsw)
        push!(solar_data, item_solar)
    end

    ## Now, we want to stitch together the data for each file
    scenario_load = reduce(vcat, load_data)
    scenario_onsw = reduce(vcat, onsw_data)
    scenario_ofsw = reduce(vcat, ofsw_data)
    scenario_solar = reduce(vcat, solar_data)

    ## Now, we want to replace this data in the files that we have
    # NOTE: ALWAYS READ IN THE 23 YEAR DATA so you can adjust the csvs accordingly
    # NOTE: Every csv EXCEPT Generators_data will be 23 years.
    # Generators_data.csv needs to be 1 yr data so costs can be multiplied

    gen_var_df = CSV.read(joinpath(base_files, "Generators_variability.csv"), DataFrame)
    load_df = CSV.read(joinpath(base_files, "Load_data.csv"), DataFrame)
    fuels_df = CSV.read(joinpath(base_files, "Fuels_data.csv"), DataFrame)
    gen_data_df = CSV.read(joinpath(base_files, "Generators_data.csv"), DataFrame)

    ## First, edit the Generators_variability.csv file
    # Only keep lengh * 8760 rows
    gen_var_df = gen_var_df[1:(scenario_length*8760), :]
    # Replace the data in the file with the scenario data
    gen_var_df[!, :solar_pv_2] = scenario_solar
    gen_var_df[!, :onshore_wind_2] = scenario_onsw
    gen_var_df[!, :fixed_offshore_wind_2] = scenario_ofsw
    gen_var_df[!, :float_offshore_wind_2] = scenario_ofsw

    ## Next, edit the Load_data.csv file
    # Only keep length * 8760 rows
    load_df = load_df[1:(scenario_length*8760), :]
    # Replace the data in the file with the scenario data
    load_df[!, :Load_MW_z2] = scenario_load
    # Replace the first value in rep_periods with the length of the scenario
    load_df[1, :Rep_Periods] = scenario_length
    # Create the Sub_Weights column (length of scenario)
    load_df[:, :Sub_Weights] .= missing
    load_df[1:scenario_length, :Sub_Weights] .= 8760


    ## Next, edit the Fuels_data.csv file
    # Only keep length * 8760 rows
    fuels_df = fuels_df[1:(scenario_length*8760), :]

    ## Next, edit the Generators_data.csv file
    # We need to multiply the costs by the length of the scenario
    gen_data_df[!, [:Inv_Cost_per_MWyr, :Inv_Cost_per_MWhyr, :Fixed_OM_Cost_per_MWyr, :Fixed_OM_Cost_per_MWhyr]] = gen_data_df[!, [:Inv_Cost_per_MWyr, :Inv_Cost_per_MWhyr, :Fixed_OM_Cost_per_MWyr, :Fixed_OM_Cost_per_MWhyr]] .* scenario_length

    ### Leave some code here to make a TDR_Results folder too
    if scenario_length != 1
        ## Grab the Period_map file
        period_map_df = CSV.read(joinpath(base_files, "TDR_Results/Period_map.csv"), DataFrame)
        # only keep length of scenario rows
        period_map_df = period_map_df[1:scenario_length, :]

        ## Make a TDR_Results folder within the scenario folder
        tdr_path = joinpath(scen_dirname_path, "TDR_Results")
        mkdir(tdr_path)

        ## Paste the Load_data, Generators_variability, Fuels_data, and Period_map files
        CSV.write(joinpath(tdr_path, "Load_data.csv"), load_df)
        CSV.write(joinpath(tdr_path, "Generators_variability.csv"), gen_var_df)
        CSV.write(joinpath(tdr_path, "Fuels_data.csv"), fuels_df)
        CSV.write(joinpath(tdr_path, "Period_map.csv"), period_map_df)
    end
    
    ### Grab any other files that you need to paste into the scenario folder
    crm_df = CSV.read(joinpath(base_files, "Capacity_reserve_margin.csv"), DataFrame)
    co2_df = CSV.read(joinpath(base_files, "CO2_cap.csv"), DataFrame)
    min_df = CSV.read(joinpath(base_files, "Minimum_capacity_requirement.csv"), DataFrame)
    net_df = CSV.read(joinpath(base_files, "Network.csv"), DataFrame)
    res_df = CSV.read(joinpath(base_files, "Reserves.csv"), DataFrame)
    
    ### Write all the necessary files to the scenario folder
    CSV.write(joinpath(scen_dirname_path, "Capacity_reserve_margin.csv"), crm_df)
    CSV.write(joinpath(scen_dirname_path, "CO2_cap.csv"), co2_df)
    CSV.write(joinpath(scen_dirname_path, "Minimum_capacity_requirement.csv"), min_df)
    CSV.write(joinpath(scen_dirname_path, "Network.csv"), net_df)
    CSV.write(joinpath(scen_dirname_path, "Reserves.csv"), res_df)
    CSV.write(joinpath(scen_dirname_path, "Generators_data.csv"), gen_data_df)
    CSV.write(joinpath(scen_dirname_path, "Generators_variability.csv"), gen_var_df)
    CSV.write(joinpath(scen_dirname_path, "Fuels_data.csv"), fuels_df)
    CSV.write(joinpath(scen_dirname_path, "Load_data.csv"), load_df)    
  
    ### Also remember to change things in the runfile for CO2_Capperiods
end