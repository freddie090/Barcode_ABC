
# Copyright 2025 Cancer Research Technology and The Institute of Cancer Research.
#
# Licensed under a software academic use license provided with this software package (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at: https://github.com/freddie090/Barcode_ABC
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

################################################################################

# Load in Euler's
e = Base.MathConstants.e

using Distributions
using DataFrames
using RCall
using CSV
using StatsBase
using Base.Threads
using Dates
using FreqTables
using DifferentialEquations

# Import Threads functions for multi-threading.
import Base.Threads.@threads
import Base.Threads.@spawn

# Create the structures and functions for creating the CancerCell
# objects: 

############
# Structures
############

# Barcode library - vector of unique barcodes to sample from when 'infecting'
# cells with barcodes.

mutable struct BarcodeLibrary
    barcodes::Array{Float64}
end

# Cancer Cell - cell_ID, barcode, birth and death rates and resistance score.

mutable struct CancerCell
    cell_ID::Int64 # ID if want to distinguish unique lineage from barcode ID.
    barcode::Float64 # Barcode identity.
    b::Float64 # birth rate.
    d::Float64 # death rate.
    R::Float64 # Resistant phenotype (binary).
    E::Float64 # Escape mutation (binary) - if 1.0, R no longer incurs cost.
    Ndiv::Int64 # Number of divisions a given cell lineage has experienced.
    #CancerCell() = new()
end

# Output of grow_cells.

mutable struct Grow_Out
    cells::Array{CancerCell}
    Nvec::Array{Int64}
    tvec::Array{Float64}
    Rvec::Array{Int64}
    Evec::Array{Int64}
    Pvec::Array{Int64}
    fin_t::Float64
    lin_df::Union{DataFrame, Nothing}
    # Constructor with default value for lin_df
    function Grow_Out(cells, Nvec, tvec, Rvec, Evec, Pvec, fin_t; lin_df=nothing)
        new(cells, Nvec, tvec, Rvec, Evec, Pvec, fin_t, lin_df)
    end
end

# Output of expand split

mutable struct Exp_Split_Out
    cells::Vector{Array{CancerCell}}
    lin_df::Union{DataFrame, Nothing}
    # Constructor with default value for lin_df
    function Exp_Split_Out(cells; lin_df=nothing)
        new(cells, lin_df)
    end
end


#############################
# Drug Concentration Function
#############################

# Given some treatment times, drug diffusion constant (k) and drug 
# strength parameter (Dc), calculate the concentration of the drug at time 
# t from t0 to tmax by solving the ODE

function drug_treat_concs(t0::Float64, tmax::Float64, k::Float64, Dc::Float64,
    treat_ons::Array{Float64}, treat_offs::Array{Float64}, dt_save_at::Float64)

    p = [Dc, 0.0]
    u0 = [0.0, 0.0]
    
    tspan = (0.0, tmax)

    function ode_fxn!(du, u, p, t)

        Dc,kp = p

        du[1] = kp
        u[2] = u[1] * Dc

    end

    prob = ODEProblem(ode_fxn!, u0, tspan, p)

    # On
    ####

    condition1(u, t, integrator) = t ∈ treat_ons
    function affect1!(integrator)
        integrator.p[2] = k 
    end

    # Off
    #####

    condition2(u, t, integrator) = t ∈ treat_offs
    function affect2!(integrator)
        integrator.p[2] = -k
    end

    # Don't let above 1.0
    condition3a(u, t, integrator) = (u[1] - 1.0)
    function affect3a!(integrator)
        integrator.u[1] = 1.0  
        integrator.p[2] = 0.0 
    end

    # Don't let below 0.0
    condition3b(u, t, integrator) = (u[1])
    function affect3b!(integrator)
        integrator.u[1] = 0.0  
        integrator.p[2] = 0.0 
    end

    cb1 = DiscreteCallback(condition1,affect1!, 
                            save_positions=(false, true))
    cb2 = DiscreteCallback(condition2,affect2!, 
                            save_positions=(false, true))
    cb3a = ContinuousCallback(condition3a,affect3a!, 
                            save_positions=(false, true))
    cb3b = ContinuousCallback(condition3b,affect3b!, 
                            save_positions=(false, true))

    cbs = CallbackSet(cb1, cb2, cb3a, cb3b)

    sol = solve(prob, callback = cbs, 
                tstops = collect(0.0:dt_save_at:tmax))

    return Dict("t"=>sol.t, 
                "dconc"=> map(x -> x[2], sol.u),
                "rconc"=> map(x -> x[1], sol.u))

end

# Return the position id of a vector, xvec, which contains the value 
# closest to the given value, x. 

function find_closest(x, xvec)
    
    n = length(xvec)
    idx = searchsortedfirst(xvec, x)
    
    if idx == 1
        return 1
    elseif idx == n
        return n
    else
        dist_prev = x - xvec[idx-1]
        dist_next = xvec[idx] - x
        return dist_prev < dist_next ? idx-1 : idx
    end
end

# Given a current time, return the drug concentration closest to that time. 

function curr_dc(curr_t::Float64, drug_concs::Dict)

    curr_dconc = drug_concs["dconc"][find_closest(curr_t, drug_concs["t"])]

    return curr_dconc

end

# ... and do the same for the relative concentration: 

function curr_rc(curr_t::Float64, drug_concs::Dict)

    curr_rconc = drug_concs["rconc"][find_closest(curr_t, drug_concs["t"])]

    return curr_rconc

end


############################
# Seed Cells Function.
############################

function seed_cells(N::Int64, b::Float64, d::Float64, rho::Float64,
    mu::Float64, sig::Float64, del::Float64, al::Float64;
    psi=0.0::Float64,
    bc_lib=BarcodeLibrary(1:N)::BarcodeLibrary, use_lib=false::Bool)

    cells = Array{CancerCell}(undef, N)
    
    # Initially, set the resistant and escape phenotypes to 0.0
    Rs = zeros(N)
    Es = zeros(N)

    # Use lib to decide whether to use a BarcodeLibrary or default to 1:N.
    if use_lib
        bcs = sample(bc_lib.barcodes, N, replace = true)
    else
        bcs = bc_lib.barcodes
    end
    # Assign all the fields to the CancerCell accordingly.
    for i in 1:length(cells)
        cells[i] = CancerCell(i, bcs[i], b, d, Rs[i], Es[i], 0)
    end

    # Assign resistance to the initial cells for the 'pre-existing' scenario
    if rho > 0.0
        nR = Int64(floor(rho*length(cells)))
        R_cells = sample(1:length(cells), nR, replace = false)
        for i in 1:nR
            cells[R_cells[i]].R = 1.0
        end
    end

    return cells

end



# Opposed to deepcopying the cell following a birth every time, it is much
# quicker to create a seperate copycell function, and call this within
# the birth-death function.

function copycell(orig_cell::CancerCell)
    new_cell::CancerCell = CancerCell(
    copy(orig_cell.cell_ID),
    copy(orig_cell.barcode),
    copy(orig_cell.b),
    copy(orig_cell.d),
    copy(orig_cell.R),
    copy(orig_cell.E),
    copy(orig_cell.Ndiv)
    )
end

###########################
# Data Collection Functions
###########################

# Replace missings in a DataFrame.

function rep_missings(df::DataFrame)

    for col in names(df)
        df[ismissing.(df[:,col]), col] .= 0.0
    end

    return df

end


# Full join an array of dataframes by colname and return, removing missings.

function join_dfs(dfs::Array{DataFrame}, colname::String)

    length(dfs) >= 2 || error("Provide >1 DataFrame in the dfs array.")

    full_df = dfs[1]
    for i in 2:length(dfs)
        full_df = outerjoin(full_df, dfs[i], on = Symbol(colname))
    end
    full_df = rep_missings(full_df)

    return full_df

end


# Barcode counts given a rep name from cells.

function get_counts(cells::Array{CancerCell}, rep_name::String)
    # return a dummy dataframe if no cells 
    if length(cells) == 0 
        df = DataFrame(bc = 0, n = 0)
    else
    df = DataFrame(bc = map(x -> x.barcode, cells))
        # depreciated by 1.4.1
        # df = by(df, :barcode, nrow)
        df = combine(groupby(df, :bc), nrow)
    end

    colnames = [:bc, Symbol(rep_name)]
    rename!(df, colnames)
    return df

end

# Barcode counts by phenotype given a rep name from cells.

function get_pheno_counts(cells::Array{CancerCell}, rep_name::String,
    top_x::Int64)

    # return a dummy dataframe if no cells 
    if length(cells) == 0 
        sub_df = DataFrame(bc = 0, nS = 0, nR = 0, nE = 0)
    else
        bcs = map(x -> x.barcode, cells)
        nEs = map(x -> x.E, cells)
        nRs = map(x -> x.R, cells)
        nSs = Float64.((nRs .== 0.0) .* (nEs .== 0.0))
        df = DataFrame(bc = map(x -> x.barcode, cells),
                       nS = nSs,
                       nR = nRs,
                       nE = nEs)
        grouped_df = combine(groupby(df, :bc), 
                             :nS => sum, 
                             :nR => sum, 
                             :nE => sum)

        grouped_df[:,:N] .= (grouped_df.nS_sum .+ grouped_df.nR_sum .+ grouped_df.nE_sum)
    end

    # Only keep the top x proportion of cells as denoted by top_prop: 
    sorted_df = sort(grouped_df, :N, rev=true)

    if nrow(sorted_df) > top_x
        sub_df = sorted_df[1:top_x ,:]
    else
        sub_df = sorted_df
    end

    sub_df[!,"rep"] .= rep_name

    colnames = vcat(:bc, :nS, :nR, :nE, :N, :rep)
    rename!(sub_df, colnames)
    return sub_df

end

# Update tracking vectors in the growth/passage function 

function update_track_vec!(kmc_out, 
    Nvec::Array{Int64}, 
    nS_vec::Array{Int64}, nR_vec::Array{Int64}, nE_vec::Array{Int64},
    tvec::Array{Float64}, Pvec::Array{Int64})

    append!(Nvec, kmc_out.Nvec)
    append!(nS_vec, kmc_out.Nvec .- kmc_out.Rvec .- kmc_out.Evec)
    append!(nR_vec, kmc_out.Rvec)
    append!(nE_vec, kmc_out.Evec)
    append!(tvec, kmc_out.tvec)
    append!(Pvec, kmc_out.Pvec)

end

# Update tracking vectors in the growth/passage function for the lineage 
# tracking version

function update_track_vec_lin_track!(kmc_out, 
    Nvec::Array{Int64}, 
    nS_vec::Array{Int64}, nR_vec::Array{Int64}, nE_vec::Array{Int64},
    tvec::Array{Float64}, Pvec::Array{Int64},
    lin_df::DataFrame)

    append!(Nvec, kmc_out.Nvec)
    append!(nS_vec, kmc_out.Nvec .- kmc_out.Rvec .- kmc_out.Evec)
    append!(nR_vec, kmc_out.Rvec)
    append!(nE_vec, kmc_out.Evec)
    append!(tvec, kmc_out.tvec)
    append!(Pvec, kmc_out.Pvec)
    append!(lin_df, kmc_out.lin_df)

end



###########################
# Growth Functions
###########################

function grow_kill_lin_kmc(cells::Array{CancerCell},
    b::Float64, d::Float64,
    mu::Float64, sig::Float64, del::Float64, al::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t0::Float64, tmax::Float64, Nmax::Int64, Cc::Int64,
    treat_ons::Array{Float64}, treat_offs::Array{Float64}, dt_save_at::Float64;
    R_real="b"::String, t_frac=0.050::Float64, treat::Bool=false,
    Passage::Int64=1, rep::String = "REP", top_x::Int64=0,
    lin_track::Bool=false)

    out_cells = deepcopy.(cells)

    0 <= mu <= 1.0 || error("mu must be between 0 and 1.")
    0 <= sig <= 1.0 || error("sig must be between 0 and 1.")
    0 <= del <= 1.0 || error("del must be between 0 and 1.")
    0 <= al <= 1.0 || error("al must be between 0 and 1.")
    0 <= t_frac <= 1.0 || error("t_frac must be between 0 and 1.")
    0 <= (al+sig) <= 1.0 || error("al and sig musn't sum to > 1.0")

    length(unique(map(x -> x.b, out_cells))) == 1 || error("For this simulation, cells should have a uniform base birth and death rate.")
    length(unique(map(x -> x.d, out_cells))) == 1 || error("For this simulation, cells should have a uniform base birth and death rate.")
    R_real == "b" || R_real == "d" || R_real == "l" || error("R_real can only take 'b', 'd' or 'l' as a value.")

    # Calculate the drug concentrations using the drug concentration ODE fxn. 
    drug_concs = drug_treat_concs(t0, tmax, k, Dc, treat_ons, treat_offs, dt_save_at)

    # Maximum cell birth and death rates - will depend on the drug 
    # concentrations - but we know what these will be based on the drug_concs 
    # vectors. 

    bmax = maximum(map(x -> x.b, out_cells))
    dmax = maximum(map(x -> x.d, out_cells))
    # Save the population's pre-cost, net growth rate, λ (lam).
    lam = bmax - dmax

    ####################################
    # Cost of Resistance Change to bdmax
    ####################################

    # If using the cost of resistance simulation, need to update bmax and dmax
    # as the highest population b/d rates can now be higher in resistant cells.
    # The cost of resistance is only implemented if δ > 0.0.
    if del > 0.0
        if R_real == "d"
            dmax = dmax + (lam*del)
            c_bdmax = bmax + dmax
        elseif R_real == "b"
            # Both are now smaller, so don't need to change either b/dmax.
            c_bdmax = bmax + dmax
        elseif R_real == "l"
            # Both are now smaller, so don't need to change either b/dmax.
            c_bdmax = bmax + dmax
        end

    ####################################
    # Non-Cost bdmax
    ####################################

    # If not using cost of resistance, simply b + d
    else
        c_bdmax = bmax + dmax
    end

    ###########################
    # Treatment Change to bdmax
    ###########################

    if treat == true

        dmax = d + (maximum(drug_concs["dconc"]))
        t_bdmax = bmax + dmax

    #####################
    # Non-Treatment bdmax
    #####################

    else
        t_bdmax = bmax + dmax
    end

    bdmax = maximum([c_bdmax, t_bdmax])
    
    # R_real determines if the cost is incurred in the birth or death rates, or
    # shared proportionately amongst both rates; "b", "d" and "l", respectively.

    ####################################

    # Create vector to hold time, and initiate at 0.
    t = t0
    tvec = Float64[t]
    # Record the population size (Nt) every tmax*t_frac.
    t_rec_change = tmax*t_frac
    t_rec = t_rec_change

    # Have a seperate counter that keeps track of the number of live cells,
    # as we can't now use length(out_cells) - includes the dead cells until the
    # end.
    Nt = length(out_cells)
    # Have a vector to save the number of live cells every t_rec.
    Nvec = Int64[Nt]
    # Also have vectors to save the number of cells that have a resistance
    # mutation (Rvec). Keep track of the number only when an event happens with
    # count vectors to avoid repeatedely using map (expensive).
    Rcount = sum(map(x -> x.R, out_cells))
    Rvec  = Int64[Rcount]
    # ... and the same for the escape mutations.
    Ecount = sum(map(x -> x.E, out_cells))
    Evec = Int64[Ecount]
    # And a vector that records the current passage of the simulation, 
    Pvec = Int64[Passage]

    if lin_track == true
        # A dataframe that will hold the top x phenotype lineage counts every 
        # t_frac:
        lin_df = DataFrame()
    end

    ##############################
    # Sampling and Killing Vectors
    ##############################

    # Opposed to sampling directly from out_cells, have a 'samp_vec'. Update
    # this following a birth or death accordingly:
    # if a birth, +1 to the length of the vector.
    # if a death, change this index in the sampling vec to 0, and add the index
    # to the 'kill_vec'.
    # Keep sampling from samp_vec until ran_cell_pos > 0
    samp_vec = collect(1:length(out_cells))
    kill_vec = Array{Int64}(undef, 0)

    ###########################
    # Grow cells until t = tmax
    ###########################

    while t <= tmax

        ran_samp_pos = rand(1:length(samp_vec))
        ran_cell_pos = samp_vec[ran_samp_pos]
        # Keep sampling from samp_vec until hit a +ve number.
        if ran_cell_pos == 0
            while ran_cell_pos == 0
                ran_samp_pos = rand(1:length(samp_vec))
                ran_cell_pos = samp_vec[ran_samp_pos]
            end
        end

        ran = rand(Uniform(0, bdmax))

        # Update the time pre-event, to prevent including events that will
        # have occured post-tmax.

        dt = -1/(bdmax * Nt) * log(rand())
        t += dt

        # If the time-change has crossed t_rec, then update the tvec and
        # Nvec tracking vectors, and update t_rec by t_rec_change.
        if t >= t_rec
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Rvec, Rcount)
            push!(Evec, Ecount)
            push!(Pvec, Passage)

            if lin_track == true
                # Get the top x phenotype counts and append to the lin_df.
            temp_lin_df = get_pheno_counts(out_cells[samp_vec[samp_vec.>0]],
                                           rep,
                                           top_x)
                temp_lin_df[:,:t] .= round(t, digits=5)
                lin_df = vcat(lin_df, temp_lin_df)
            end
            t_rec += t_rec_change
        end

        # If the time-change has crossed tmax, then there wasn't time for
        # the birth or event to happen, t is capped at tmax, and the
        # simulation ends.
        if t > tmax
            t = tmax
            break
        end
        # If N is now >= Nmax, also end the simulation, whilst updating the
        # tracking vectors.
        if Nt >= Nmax
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Rvec, Rcount)
            push!(Evec, Ecount)
            push!(Pvec, Passage)
            break
        end

        # Otherwise, continue with the birth and death steps.

        #######################################
        # Realisation of the Cost of Resistance
        #######################################

        # NB - the cost of resistance is realied as follows,
        # where λ = (b - d).

        # trade-off realised in birth rate.
            # b = b - λ*δ
            # d = d
        # trade-off realised in death rate.
            # b = b
            # d = d + λ*δ
        # trade-off realised in both rates.
            # b = b * (1-δ)
            # d = d * (1-δ)

        # Set cell_b and cell_d using the trade-off if performing the cost of
        # resistance simulation.

        if del > 0.0

            # If chosen cell is resistant...
            if out_cells[ran_cell_pos].R > 0

                # Calculate the cost incurred according to 'R_real'.
                if R_real == "b"
                    cell_b = out_cells[ran_cell_pos].b - (lam*del)
                    cell_d = out_cells[ran_cell_pos].d
                end
                if R_real == "d"
                    cell_b = out_cells[ran_cell_pos].b
                    cell_d = out_cells[ran_cell_pos].d + (lam*del)
                end
                if R_real == "l"
                    cell_b = out_cells[ran_cell_pos].b * (1-del)
                    cell_d = out_cells[ran_cell_pos].d * (1-del)
                end

            # Otherwise, cell is sensitive or escape...
            else
                cell_b = out_cells[ran_cell_pos].b
                cell_d = out_cells[ran_cell_pos].d
            end

            # Set to 0 if < 0
            if cell_b < 0
                cell_b = 0
            end
            if cell_d < 0
                cell_d = 0
            end

        #######################################

        # If not using cost of resistance simulation, simply assign rates from
        # the randomly selected cell.
        else
            cell_b = out_cells[ran_cell_pos].b
            cell_d = out_cells[ran_cell_pos].d
        end

        #########################
        # Treatment Induced Death
        #########################

        if treat == true

            # Get current drug concentration:
            curr_dconc = curr_dc(t, drug_concs)

            # If chosen cell is resistant:
            if out_cells[ran_cell_pos].R > 0

                cell_d = cell_d + (curr_dconc * (1 - psi))

            # If chosen cell is escape:
            elseif out_cells[ran_cell_pos].E > 0

                cell_d = cell_d + (curr_dconc * (1 - psi))

            # Otherwise, cell is sensitive...
            else

                cell_d = cell_d + (curr_dconc)

            end

        end

        #################
        # Logistic Growth
        #################

        # Modify growth and death rates according to carrying capacity (Cc) - 
        # the reduction in net growth rate is distributed evently between the 
        # birth and death rate. 

        cell_b = cell_b * (1 - (Nt/Cc))
        cell_d = cell_d * (1 - (Nt/Cc))

        #############################
        # Birth step
        #############################

        if ran < cell_b

            # Update the cells number of divisions.
            out_cells[ran_cell_pos].Ndiv += 1
            ran_cell = copycell(out_cells[ran_cell_pos])

            (ran_cell.R + ran_cell.E) <= 1.0 || error("A cell's R and E phenotype scores should not sum to > 1.0.")

            ran_p = rand()

            ########
            # Escape
            ########
            
            # just birth
            if ran_cell.E > 0.0
                Ecount += 1

            #############################
            # Sensitive -> Resistant
            #############################

            elseif ran_cell.R == 0.0
                if mu > ran_p
                        ran_cell.R += 1.0
                        # Also update Rcount with the new number of R cells.
                        Rcount += 1
                end

            ##################################
            # Resistant -> Sensitive or Escape
            ##################################

            elseif ran_cell.R > 0.0

                # Get current relative drug concentration to make escape 
                # transitions drug-dependent:
                curr_rconc = curr_rc(t, drug_concs)
                # Scale alpha by relative concentration: 
                al_scal = al*curr_rconc

                # R -> S
                if ran_p < sig 
                    ran_cell.R = 0.0
                # R -> E
                elseif sig <= ran_p < (sig + al_scal)
                    ran_cell.R -= 1.0
                    ran_cell.E += 1.0
                    Ecount += 1
                # just birth
                else 
                    Rcount += 1
                end
            
            end

            #############################
            # Update count vectors
            #############################

            # Add the updated daughter cells to the cell vector. 
            push!(out_cells, ran_cell)
            # Update samp_vec with the new length of the output cells.
            push!(samp_vec, (length(out_cells)))
            # Update number of cells to normalise dt.
            Nt += 1

        end

        #############################
        # Death step
        #############################

        if bmax <= ran < bmax + cell_d
            # Remove this chosen cell's index from the samp_vec.
            samp_vec[ran_samp_pos] = 0
            push!(kill_vec, ran_cell_pos)
            # Update number of cells to normalise dt.
            Nt -= 1
            # Update the R and E count vectors accordingly.
            if out_cells[ran_cell_pos].R == 1.0
                Rcount -= 1
            elseif out_cells[ran_cell_pos].E == 1.0
                Ecount -= 1
            end
        end

        #############################
        # Check for extinction
        #############################

        # Break if no cells remain. Make use of the fact that:
        # samp_vec = N0 + nbirths
        # kill_vec = ndeaths
        # therefore if length(kill_vec) >= length(samp_vec), there are no cells
        # remaining. Also save this outcome to the record vectors, removing 
        # the previous values from the record vector to avoid multiple times.

        if length(kill_vec) >= length(samp_vec)

            deleteat!(Nvec, length(Nvec))
            deleteat!(tvec, length(tvec))
            deleteat!(Rvec, length(Rvec))
            deleteat!(Evec, length(Evec))
            deleteat!(Pvec, length(Pvec))

            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Rvec, Rcount)
            push!(Evec, Ecount)
            push!(Pvec, Passage)

            break
        end

    end

    # Now perform all the kill steps using the indexes stored in kill_vec.
    deleteat!(out_cells, sort(kill_vec))

    fin_t = round(t, digits = 4)

    if lin_track == true
        return Grow_Out(out_cells, Nvec, tvec, Rvec, Evec, Pvec, fin_t; 
                        lin_df=lin_df)
    else
        return Grow_Out(out_cells, Nvec, tvec, Rvec, Evec, Pvec, fin_t)
    end
    
end


##################################
# Expand, Split Cells Function.
##################################

# Expand N cells, split equally into  N_reps replicates, each containing
# N_seed cells. 

function expand_split_cells(N::Int64, b::Float64, d::Float64,
    rho::Float64, mu::Float64, sig::Float64, del::Float64, al::Float64,
    t_exp::Float64, N_seed::Int64, N_reps::Int64;
    psi=0.0::Float64, t_frac=0.050::Float64,
    R_real="b"::String, 
    top_x::Int64 = 0, lin_track::Bool=false)

    # Use limiting probabilities according to use_lim_probs, and plastic/res_cost
    init_cells = seed_cells(N, b, d, rho, mu, sig, del, al, psi=psi)

    # Expand the cells for t_exp.
    exp_out = grow_kill_lin_kmc(init_cells, b, d, mu, sig, del, al, 
                                0.0, 0.0, 0.0,
                                0.0, t_exp, 10^10, 10^10,
                                [0.0], [0.0], 1e-01,
                                t_frac=t_frac, rep="POT",
                                R_real=R_real, treat=false,
                                top_x=top_x,
                                lin_track=lin_track)

    exp_cells = exp_out.cells
    if lin_track == true
        exp_lin_df = exp_out.lin_df
    end

    # Need to split into N_reps
    N_reps * N_seed < length(exp_cells) || error("There are not enough cells to split into the chosen number of replicates.")

    # Sample all the cells at once, and then reshape into replicates.
    rep_cells = sample(exp_cells, (N_reps * N_seed), replace = false)
    rep_cells = reshape(rep_cells, (N_seed, N_reps))

    # Return as a vector of cell objects: 
    fin_rep_cells = Array{Vector{CancerCell}}(undef, N_reps)
    for i in 1:N_reps
        fin_rep_cells[i] = rep_cells[:,i]
    end

    if lin_track == true
        return Exp_Split_Out(fin_rep_cells, lin_df=exp_lin_df)
    else
        return Exp_Split_Out(fin_rep_cells)
    end

end


#########################
# Drug Kill Passage Cells
#########################

# A function that grows and kills the cells whilst also passaging them either 
# when t >= t_Pass, or N >= Nmax.

function grow_kill_pass_lin_kmc(cells::Array{CancerCell}, 
    n0::Int64, b::Float64, d::Float64, 
    mu::Float64, sig::Float64, del::Float64, al::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t0::Float64, tmax::Float64, t_Pass::Float64, Nmax::Int64, Cc::Int64,
    treat_ons::Array{Float64}, treat_offs::Array{Float64}, dt_save_at::Float64;
    R_real="b"::String, t_frac=0.005::Float64, treat::Bool=false,
    n_Pass::Int64=2, rep::Int64=0,
    top_x::Int64=0, lin_track::Bool=false)

    # Set initial time, Passage time and Passage id
    
    curr_t = t0

    if tmax < t_Pass
        next_t = tmax
    else
        next_t = t_Pass
    end

    curr_P = 1

    # Vector to hold the cell lineage dataframes. 
    cell_lin_df_vec = Array{DataFrame}(undef, 0)

    if lin_track == true
        # A dataframe that holds the top x barcodes every t_rec: 
        sub_lin_df = DataFrame()
    end

    # Vectors to hold the total population (N), time (t), Passage (P), 
    # sensitive (nS) and resistant (nR) population vectors. 

    Nvec = Int64[]
    nS_vec = Int64[]
    nR_vec = Int64[]
    nE_vec = Int64[]
    tvec = Float64[]
    Pvec = Int64[]
    
    while curr_t <= tmax

        # Grow/kill cells: 

        kmc_out = grow_kill_lin_kmc(cells, b, d, mu, sig, del, al,
                                    Dc, k, psi, curr_t, next_t, Nmax, Cc,
                                    treat_ons, treat_offs, dt_save_at,
                                    R_real=R_real, t_frac=t_frac,
                                    treat=treat,Passage=curr_P,
                                    rep = string("DT", rep), top_x=top_x,
                                    lin_track=lin_track)
                                    
        # Check for Nmax exceeded 
        #########################
        if last(kmc_out.Nvec) >= Nmax
            # If already on n_Pass, end simulation here
            if last(kmc_out.Pvec) >= n_Pass
                # Update cell linege dataframe vector with these cells
                bc_df = get_counts(kmc_out.cells, string("DT", rep, "_P", curr_P))
                push!(cell_lin_df_vec, bc_df)
                # Update tracker vectors
                if lin_track == true
                    update_track_vec_lin_track!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec, sub_lin_df)
                else
                    update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
                end
                # Update current time
                curr_t = last(kmc_out.tvec)
                # ... and terminate
                println("Population exceeded Nmax.")
                break
            # Otherwise continue with Passaging cells
            else
                # Update cell linege dataframe vector with these cells
                bc_df = get_counts(kmc_out.cells, string("DT", rep, "_P", curr_P))
                push!(cell_lin_df_vec, bc_df)
                # Update tracker vectors
                if lin_track == true
                    update_track_vec_lin_track!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec, sub_lin_df)
                else
                    update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
                end
                # Check enough cells to sample
                if length(kmc_out.cells) < n0
                    println("Not enough cells to sample from.")
                    break
                else
                    # Sample from these cells to Passage
                    cells = sample(kmc_out.cells, n0, replace=false)
                    # Update Passage number
                    curr_P += 1
                    # Update time and continue with simulation 
                    curr_t = last(kmc_out.tvec)
                    next_t = tmax                    
                    # Continue with simulation 
                end
            end

        # Check for population extinction 
        #################################
        elseif length(kmc_out.cells) == 0
            # Update cell linege dataframe vector with these cells
            bc_df = get_counts(kmc_out.cells, string("DT", rep, "_P", curr_P))
            push!(cell_lin_df_vec, bc_df)
            # Update tracker vectors
            if lin_track == true
                update_track_vec_lin_track!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec, sub_lin_df)
            else
                update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
            end
            # Update current time
            curr_t = last(kmc_out.tvec)
            # ... and terminate
            println("Population reached extinction.")
            break

        # Check for crossing tmax
        #########################
        elseif last(kmc_out.tvec) >= tmax
            # Update cell linege dataframe vector with these cells
            bc_df = get_counts(kmc_out.cells, string("DT", rep, "_P", curr_P))
            push!(cell_lin_df_vec, bc_df)
            # Update tracker vectors
            if lin_track == true
                update_track_vec_lin_track!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec, sub_lin_df)
            else
                update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
            end
            # Update current time
            curr_t = last(kmc_out.tvec)
            # ... and terminate
            println("tmax exceeded.")
            break

        # Check for crossing Passage time
        #################################
        elseif last(kmc_out.tvec) >= t_Pass
            # Only perform Passaging if not already on n_Pass
            if last(kmc_out.Pvec) < n_Pass
                # Check enough cells to sample
                if length(kmc_out.cells) < n0
                    println("Not enough cells to sample from.")
                    # Update cells
                    cells = kmc_out.cells
                    # Update tracker vectors
                    if lin_track == true
                        update_track_vec_lin_track!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec, sub_lin_df)
                    else
                        update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
                    end
                        # Update time and continue with simulation 
                    curr_t = last(kmc_out.tvec)
                    next_t = tmax
                else
                    # Update cell linege dataframe vector with these cells
                    bc_df = get_counts(kmc_out.cells, string("DT", rep, "_P", curr_P))
                    push!(cell_lin_df_vec, bc_df)
                    # Update tracker vectors
                    if lin_track == true
                        update_track_vec_lin_track!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec, sub_lin_df)
                    else
                        update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
                    end
                    # Sample from these cells to Passage
                    cells = sample(kmc_out.cells, n0, replace=false)
                    # Update Passage number
                    curr_P += 1
                    # Update time and continue with simulation 
                    curr_t = last(kmc_out.tvec)
                    next_t = tmax
                end
            end
        end
    end

    # If not all Passages realised, repeat the most recent Passage which was.
    for i in 1:n_Pass
        if i != 1
            if length(cell_lin_df_vec) < i
                temp_df = deepcopy(cell_lin_df_vec[i-1])
                colnames = [:bc, Symbol("DT", rep, "_P", i)]
                rename!(temp_df, colnames)
                push!(cell_lin_df_vec, temp_df)
            end
        end
    end
    
    if lin_track == true
        return Dict("cell_lin_df_vec"=>cell_lin_df_vec, 
                    "Nvec"=>Nvec, "tvec"=>tvec, "Pvec"=>Pvec, 
                    "nS_vec"=>nS_vec, "nR_vec"=>nR_vec, "nE_vec"=>nE_vec,
                    "sub_lin_df"=>sub_lin_df)
    else
        return Dict("cell_lin_df_vec"=>cell_lin_df_vec, 
                    "Nvec"=>Nvec, "tvec"=>tvec, "Pvec"=>Pvec, 
                    "nS_vec"=>nS_vec, "nR_vec"=>nR_vec, "nE_vec"=>nE_vec)
    end
end

################################################################################ 

# A version that returns the full simulation (t & N vecs) output, alongside 
# the lineage and solution dataframes. The sub-sampled lineage dataframes at
# multiple timepoints are returned if 'lin_track=true'. 

function exp_pulse_treat_kmc_full_sim_fxn(n0::Int64, b::Float64, d::Float64,
    rho::Float64, mu::Float64, sig::Float64, del::Float64, al::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t0::Float64, t_exp::Float64, tmax::Float64, t_Pass::Float64, 
    Nmax::Int64, Cc::Int64,
    treat_ons::Array{Float64}, treat_offs::Array{Float64}, dt_save_at::Float64,
    t_keep::Array{Float64};
    n_Pass::Int64 = 2,
    R_real::String="b", t_frac::Float64=0.005,
    top_x::Int64=0, lin_track::Bool=false)

    # Expand and split into 4 replicates
    rep_out = expand_split_cells(n0, b, d, rho, mu, sig, del, 0.0,
                                   t_exp, n0, 4, R_real=R_real,
                                   top_x=top_x, lin_track=lin_track)

    rep_cells = rep_out.cells

    # Vectors to keep final outputs: 
    fin_t_outs = Float64[]
    fin_u_outs = Int64[]
    df_outs = DataFrame[]
    lin_dfs = DataFrame[]
    # Vector to keep solution output dataframes
    sim_dfs = DataFrame[]

    if lin_track == true
        # Vector to keep sub-top-x lineage dataframes
        sub_lin_dfs = DataFrame[]
        # Add the POT lineages
        push!(sub_lin_dfs, rep_out.lin_df)
    end

    # Repeat simulation for each replicate
    #Threads.@threads for i in 1:4
    for i in 1:4

        # Run to derive simulation output 
        sim = grow_kill_pass_lin_kmc(rep_cells[i],
                                     n0, b, d, mu, sig, del, al,
                                     Dc, k, psi, 
                                     t0, tmax, t_Pass, Nmax, Cc,
                                     treat_ons, treat_offs, dt_save_at,
                                     R_real=R_real, t_frac=t_frac,
                                     n_Pass=n_Pass, rep=i,
                                     treat=true,
                                     top_x=top_x, lin_track=lin_track)

        # Create a dataframe containing the solution output: 
        sim_df = DataFrame(
            t = sim["tvec"], 
            N = sim["Nvec"],
            nS = sim["nS_vec"], 
            nR = sim["nR_vec"], 
            nE = sim["nE_vec"],
            rep = i
        )

        push!(sim_dfs, sim_df)

        # Combine the lineage matrices into a single dataframe
        if n_Pass > 1
            out_df = join_dfs(sim["cell_lin_df_vec"], "bc")
        else
            out_df = sim["cell_lin_df_vec"][1]
        end

        push!(df_outs, out_df)

        if lin_track==true
            # Save the replicate's top x subset lineage dataframe:
            push!(sub_lin_dfs, sim["sub_lin_df"])
        end

        # Vectors to save output t and population sizes (u) - size of 
        # t_keep vector, plus number of passages. 
        t_outs = Vector{Float64}(undef, length(t_keep)+n_Pass)
        u_outs = Vector{Int64}(undef, length(t_keep)+n_Pass)

        # Save the t_keep times. 
        for j in 1:length(t_keep)
            # If the current time wasn't in the solution, find closest 
            # value to keep. 
            if !(t_keep[j] in sim["tvec"])
                # Find the closest time reached for this value in t_keep
                t_closest_pos = findmin(abs.(sim["tvec"] .- t_keep[j]))[2]
                t_closest = sim["tvec"][t_closest_pos]
                t_outs[j] = t_closest
                # Keep the pop size reached at this time
                u_outs[j] = sim["Nvec"][t_closest_pos]
            # Otherwise keep the closest solution 
            else
                t_realised_pos = findlast(sol["tvec"] .== t_keep[j][1])
                t_realised = sol["tvec"][t_realised_pos]
                t_outs[j] = t_realised
                u_outs[j] = sol["Nvec"][t_realised_pos]
            end
        end

        # Extract the Passage times using the given integrator variable. 
        # If Passage not reached, use the most recent time. 
        for j in 1:n_Pass
            if sum(sim["Pvec"] .== j) > 0
                t_realised_pos = findlast(sim["Pvec"] .== j)
                t_outs[length(t_keep)+j] = sim["tvec"][t_realised_pos]
                u_outs[length(t_keep)+j] = sim["Nvec"][t_realised_pos]
            else
                # If current passage not realised, save previous pop size. 
                t_outs[length(t_keep)+j] = t_outs[length(t_keep)+j-1]
                u_outs[length(t_keep)+j] = u_outs[length(t_keep)+j-1]
            end
        end

        append!(fin_t_outs, t_outs)
        append!(fin_u_outs, u_outs)

    end

    # Round t and pop sizes
    fin_t_outs = Float64.(round.(fin_t_outs, digits = 0))

    # Combine the lineage DataFrames
    fin_lin_df = join_dfs(df_outs, "bc")
    # Combine the solution DataFrames
    fin_sol_df = vcat(sim_dfs...)

    if lin_track == true
        return Dict("t"=>fin_t_outs, "u"=>fin_u_outs, 
                    "lin_df"=>fin_lin_df, "sol_df"=>fin_sol_df,
                    "sub_lin_dfs"=>sub_lin_dfs)
    else
        return Dict("t"=>fin_t_outs, "u"=>fin_u_outs, 
                    "lin_df"=>fin_lin_df, "sol_df"=>fin_sol_df)
    end
    
end

################################################################################



