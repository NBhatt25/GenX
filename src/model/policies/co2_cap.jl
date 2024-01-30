"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""

@doc raw"""
	co2_cap!(EP::Model, inputs::Dict, setup::Dict)

This policy constraints mimics the CO$_2$ emissions cap and permit trading systems, allowing for emissions trading across each zone for which the cap applies. The constraint $p \in \mathcal{P}^{CO_2}$ can be flexibly defined for mass-based or rate-based emission limits for one or more model zones, where zones can trade CO$_2$ emissions permits and earn revenue based on their CO$_2$ allowance. Note that if the model is fully linear (e.g. no unit commitment or linearized unit commitment), the dual variable of the emissions constraints can be interpreted as the marginal CO$_2$ price per tonne associated with the emissions target. Alternatively, for integer model formulations, the marginal CO$_2$ price can be obtained after solving the model with fixed integer/binary variables.

The CO$_2$ emissions limit can be defined in one of the following ways: a) a mass-based limit defined in terms of annual CO$_2$ emissions budget (in million tonnes of CO2), b) a load-side rate-based limit defined in terms of tonnes CO$_2$ per MWh of demand and c) a generation-side rate-based limit defined in terms of tonnes CO$_2$ per MWh of generation.

**Mass-based emissions constraint**

Mass-based emission limits are implemented in the following expression. For each constraint, $p \in \mathcal{P}^{CO_2}_{mass}$, we define a set of zones $z \in \mathcal{Z}^{CO_2}_{p,mass}$ that can trade CO$_2$ allowance. Input data for each constraint  $p \in \mathcal{P}^{CO_2}_{mass}$ requires the CO$_2$ allowance/ budget for each model zone, $\epsilon^{CO_{2}}_{z,p, mass}$, to be provided in terms of million metric tonnes. For every generator $y$, the parameter $\epsilon_{y,z}^{CO_2}$ reflects the specific $CO_2$ emission intensity in tCO$_2$/MWh associated with its operation.  The resulting constraint is given as:

```math
\begin{aligned}
    \sum_{z \in \mathcal{Z}^{CO_2}_{p,mass}} \sum_{y \in \mathcal{G}} \sum_{t \in \mathcal{T}} \left(\epsilon_{y,z}^{CO_2} \times \omega_{t} \times \Theta_{y,z,t} \right)
   & \leq \sum_{z \in \mathcal{Z}^{CO_2}_{p,mass}} \epsilon^{CO_{2}}_{z,p, mass} \hspace{1 cm}  \forall p \in \mathcal{P}^{CO_2}_{mass}
\end{aligned}
```

In the above constraint, we include both power discharge and charge term for each resource to account for the potential for CO$_2$ emissions (or removal when considering negative emissions technologies) associated with each step. Note that if a limit is applied to each zone separately, then the set $\mathcal{Z}^{CO_2}_{p,mass}$ will contain only one zone with no possibility of trading. If a system-wide emission limit constraint is applied, then $\mathcal{Z}^{CO_2}_{p,mass}$ will be equivalent to a set of all zones.

**Load-side rate-based emissions constraint**

We modify the right hand side of the above mass-based constraint, $p \in \mathcal{P}^{CO_2}_{load}$, to set emissions target based on a CO$_2$ emission rate limit in tCO$_2$/MWh $\times$ the total demand served in each zone. In the following constraint, total demand served takes into account non-served energy and storage related losses. Here, $\epsilon_{z,p,load}^{maxCO_2}$ denotes the emission limit in terms on tCO$_2$/MWh.

```math
\begin{aligned}
    \sum_{z \in \mathcal{Z}^{CO_2}_{p,load}} \sum_{y \in \mathcal{G}} \sum_{t \in \mathcal{T}} \left(\epsilon_{y,z}^{CO_2} \times \omega_{t} \times \Theta_{y,t,z} \right)
    \leq & \sum_{z \in \mathcal{Z}^{CO_2}_{p,load}} \sum_{t \in \mathcal{T}}  \left(\epsilon_{z,p,load}^{CO_2} \times  \omega_{t} \times D_{z,t} \right) \\  + & \sum_{z \in \mathcal{Z}^{CO_2}_{p,load}} \sum_{y \in \mathcal{O}}  \sum_{t \in \mathcal{T}} \left(\epsilon_{z,p,load}^{CO_2} \times \omega_{t} \times \left(\Pi_{y,t,z} - \Theta_{y,t,z} \right) \right) \\  - & \sum_{z \in \mathcal{Z}^{CO_2}_{p,load}} \sum_{s \in \mathcal{S} } \sum_{t \in \mathcal{T}}  \left(\epsilon_{z,p,load}^{CO_2} \times \omega_{t} \times \Lambda_{s,z,t}\right) \hspace{1 cm}  \forall p \in \mathcal{P}^{CO_2}_{load}
\end{aligned}
```

**Generator-side emissions rate-based constraint**

Similarly, a generation based emission constraint is defined by setting the emission limit based on the total generation times the carbon emission rate limit in tCO$_2$/MWh of the region. The resulting constraint is given as:

```math
\begin{aligned}
\sum_{z \in \mathcal{Z}^{CO_2}_{p,gen}} \sum_{y \in \mathcal{G}} \sum_{t \in \mathcal{T}} & \left(\epsilon_{y,z}^{CO_2} \times \omega_{t} \times \Theta_{y,t,z} \right) \\
    \leq \sum_{z \in \mathcal{Z}^{CO_2}_{p,gen}} \sum_{y \in \mathcal{G}} \sum_{t \in \mathcal{T}} & \left(\epsilon_{z,p,gen}^{CO_2} \times  \omega_{t} \times \Theta_{y,t,z} \right)  \hspace{1 cm}  \forall p \in \mathcal{P}^{CO_2}_{gen}
\end{aligned}
```

Note that the generator-side rate-based constraint can be used to represent a fee-rebate (``feebate'') system: the dirty generators that emit above the bar ($\epsilon_{z,p,gen}^{maxCO_2}$) have to buy emission allowances from the emission regulator in the region $z$ where they are located; in the same vein, the clean generators get rebates from the emission regulator at an emission allowance price being the dual variable of the emissions rate constraint.
"""
function co2_cap!(EP::Model, inputs::Dict, setup::Dict)

	println("C02 Policies Module")

	dfGen = inputs["dfGen"]
	SEG = inputs["SEG"]  # Number of lines
	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

	### Constraints ###

	## Mass-based: Emissions constraint in absolute emissions limit (tons)
	if setup["CO2Cap"] == 1
		if setup["CO2CapPeriods"] <= 1
			# CO2 cap based on a single period
			println("-- Using a single period for mass-based CO2 cap")
			EP[:cCO2Emissions_systemwide] = single_mass_cap(EP, inputs, 1:T)
		else
			println("-- Using $(setup["CO2CapPeriods"]) periods for mass-based CO2 cap")
			# Breaking up into the timeseries into CO2CapPeriods periods
			# If periods are uneven, make the first period longest
			period_times = timeseries2periods(T, setup["CO2CapPeriods"])
			time_ranges = [p[1]:p[2] for p in eachrow(period_times)]
			if "CO2CapWiggle" in keys(setup) && setup["CO2CapWiggle"] > 0
				# Each period can deviate by CO2CapWiggle % from the average
				println("-- Allowing $(100*setup["CO2CapWiggle"])% deviation from average CO2 cap in each period")
				EP[:cCO2Emissions_systemwide] = multi_mass_cap(EP, inputs, time_ranges, setup["CO2CapWiggle"])
				# The total must still be less than the cap
				EP[:cCO2Emissions_systemwide_tot] = single_mass_cap(EP, inputs, 1:T)
			else
			# Each period must strictly stay within the per-period limit
				EP[:cCO2Emissions_systemwide] = multi_mass_cap(EP, inputs, time_ranges)
			end
		end

	## Load + Rate-based: Emissions constraint in terms of rate (tons/MWh)
	elseif setup["CO2Cap"] == 2 ##This part moved to non_served_energy.jl
		if setup["CO2CapPeriods"] <= 1
			println("-- Using a single period for load-based emissions intensity CO2 cap")
			EP[:cCO2Emissions_systemwide] = single_intensity_cap_load(EP, inputs, 1:T, 0.0, setup["StorageLosses"])
		else
			println("-- Using $(setup["CO2CapPeriods"]) periods for load-based emissions intensity CO2 cap")
			period_times = timeseries2periods(T, setup["CO2CapPeriods"])
			time_ranges = [p[1]:p[2] for p in eachrow(period_times)]
			if "CO2CapWiggle" in keys(setup) && setup["CO2CapWiggle"] > 0
				println("-- Allowing $(100*setup["CO2CapWiggle"])% deviation from average CO2 cap in each period")
				EP[:cCO2Emissions_systemwide] = multi_intensity_cap_load(EP, inputs, time_ranges, setup["CO2CapWiggle"], setup["StorageLosses"])
				EP[:cCO2Emissions_systemwide_tot] = single_intensity_cap_load(EP, inputs, 1:T, 0.0, setup["StorageLosses"])
			else
				EP[:cCO2Emissions_systemwide] = multi_intensity_cap_load(EP, inputs, time_ranges, 0.0, setup["StorageLosses"])
			end
		end
	## Generation + Rate-based: Emissions constraint in terms of rate (tons/MWh)
	elseif (setup["CO2Cap"]==3)
		if setup["CO2CapPeriods"] <= 1

			EP[:cCO2Emissions_systemwide] = single_intensity_cap_gen(EP, inputs, 1:T)
		else
			period_times = timeseries2periods(T, setup["CO2CapPeriods"])
			time_ranges = [p[1]:p[2] for p in eachrow(period_times)]
			if "CO2CapWiggle" in keys(setup) && setup["CO2CapWiggle"] > 0
				EP[:cCO2Emissions_systemwide] = multi_intensity_cap_gen(EP, inputs, time_ranges, setup["CO2CapWiggle"])
				EP[:cCO2Emissions_systemwide_tot] = single_intensity_cap_gen(EP, inputs, 1:T)
			else
				EP[:cCO2Emissions_systemwide] = multi_intensity_cap_gen(EP, inputs, time_ranges)
			end
		end
	end 

end

function timeseries2periods(T::Int64, numperiods::Int64)
	# Breaks up the timeseries into numperiods periods
	# If periods are uneven, make the first period longest
	period_times = Array{Int64}(undef, numperiods, 2)
	min_period_length = Int64(floor(T/numperiods))
	remainder = T % numperiods
	for p in 1:numperiods
		period_times[p,1] = (p-1)*min_period_length + 1 + remainder
		period_times[p,2] = p*min_period_length + remainder
	end
	period_times[1,1] = 1
	return period_times
end

function single_mass_cap(EP::Model, inputs::Dict, time_range::AbstractRange{Int64}, period_deviation::Float64=0.0)
	return @constraint(EP, [cap=1:inputs["NCO2Cap"]],
		sum(inputs["omega"][t] * EP[:eEmissionsByZone][z,t] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]), t in time_range) 
		<=
		(1 + period_deviation) * sum(inputs["dfMaxCO2"][z,cap] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]))
	)
end

function multi_mass_cap(EP::Model, inputs::Dict, time_ranges::Vector{<:AbstractRange{Int64}}, period_deviation::Float64=0.0)
	num_periods = size(time_ranges,1)
	return @constraint(EP, [cap=1:inputs["NCO2Cap"], p=1:num_periods],
		sum(inputs["omega"][t] * EP[:eEmissionsByZone][z,t] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]), t in time_ranges[p]) 
		<=
		(1 + period_deviation) * sum(inputs["dfMaxCO2"][z,cap] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap])) / num_periods
	)
end

function single_intensity_cap_load(EP::Model, inputs::Dict, time_range::AbstractRange{Int64}, period_deviation::Float64=0.0, include_stor_losses::Int64=0)
	load_emissions = @expression(EP, [cap=1:inputs["NCO2Cap"]], 
		(1 + period_deviation) * sum(
			+ inputs["dfMaxCO2Rate"][z,cap] * sum(
				inputs["omega"][t] * (inputs["pD"][t,z] - sum(EP[:vNSE][s,t,z] for s in 1:inputs["SEG"])) 
				for t in time_range
			) 
			for z = findall(x->x==1, inputs["dfCO2CapZones"][:,cap])
		) 
	)
	if include_stor_losses == 1
		println("-- Including storage losses in load-based emissions constraint")
		storage_losses = @expression(EP, [cap=1:inputs["NCO2Cap"]], 
			(1 + period_deviation) * sum(inputs["dfMaxCO2Rate"][z,cap] * EP[:eELOSSByZone][z] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]))
		)
		add_to_expression!(load_emissions, storage_losses)
	end
	return @constraint(EP, [cap=1:inputs["NCO2Cap"]],
		sum(inputs["omega"][t] * EP[:eEmissionsByZone][z,t] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]), t in time_range) 
		<=
		load_emissions[cap]
	)
end

function multi_intensity_cap_load(EP::Model, inputs::Dict, time_ranges::Vector{<:AbstractRange{Int64}}, period_deviation::Float64=0.0, include_stor_losses::Int64=0)
	num_periods = size(time_ranges,1)
	load_emissions = @expression(EP, [cap=1:inputs["NCO2Cap"], p=1:num_periods], 
		(1 + period_deviation) * sum(
			+ inputs["dfMaxCO2Rate"][z,cap] * sum(
				inputs["omega"][t] * (inputs["pD"][t,z] - sum(EP[:vNSE][s,t,z] for s in 1:inputs["SEG"])) 
				for t in time_ranges[p]
			) 
			for z = findall(x->x==1, inputs["dfCO2CapZones"][:,cap])
		) 
	)
	if include_stor_losses == 1
		println("-- Including storage losses in load-based emissions constraint")
		storage_losses = @expression(EP, [cap=1:inputs["NCO2Cap"], p=1:num_periods], 
			(1 + period_deviation) * sum(inputs["dfMaxCO2Rate"][z,cap] * EP[:eELOSSByZone][z] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]))
		)
		add_to_expression!(load_emissions, storage_losses)
	end
	return @constraint(EP, [cap=1:inputs["NCO2Cap"], p=1:num_periods],
		sum(inputs["omega"][t] * EP[:eEmissionsByZone][z,t] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]), t in time_ranges[p])
		<=
		load_emissions[cap, p]
	)
end

function single_intensity_cap_gen(EP::Model, inputs::Dict, time_range::AbstractRange{Int64}, period_deviation::Float64=0.0)
	return @constraint(EP, [cap=1:inputs["NCO2Cap"]],
		sum(inputs["omega"][t] * EP[:eEmissionsByZone][z,t] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]), t in time_range) 
		<=
		(1 + period_deviation) * sum(inputs["dfMaxCO2Rate"][z,cap] * inputs["omega"][t] * EP[:eGenerationByZone][z,t] for t in time_range, z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]))
	)
end

function multi_intensity_cap_gen(EP::Model, inputs::Dict, time_ranges::Vector{<:AbstractRange{Int64}}, period_deviation::Float64=0.0)
	num_periods = size(time_ranges,1)
	return @constraint(EP, [cap=1:inputs["NCO2Cap"], p=1:num_periods],
		sum(inputs["omega"][t] * EP[:eEmissionsByZone][z,t] for z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]), t in time_ranges[p]) 
		<=
		(1 + period_deviation) * sum(inputs["dfMaxCO2Rate"][z,cap] * inputs["omega"][t] * EP[:eGenerationByZone][z,t] for t in time_ranges[p], z=findall(x->x==1, inputs["dfCO2CapZones"][:,cap]))
	)
end