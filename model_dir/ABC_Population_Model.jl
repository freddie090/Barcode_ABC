
# Copyright 2025 Cancer Research Technology and The Institute of Cancer Research.
#
# Licensed under a software academic use license provided with this software package (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at: https://github.com/freddie090/Barcode_ABC
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

################################################################################

module ABC_Population_Model

# This module is used within the pyABC workflow to fit a model of 
# evolution of 2/3 (depending on model A,B,C) phenotypic compartments through
# treatment to observed data - in this first step, this consists of total
# cell size populations at given timepoints.
# The model combines a deterministic ODE with stochastic jump process, switching
# between the two when a given threshold is crossed.

using DifferentialEquations
using JumpProcesses
using Distributions
using Suppressor
using RCall
using Distances
using Hungarian
using Combinatorics

e = Base.MathConstants.e

pop_fun(x) = return(x[2:4])
pass_fun(x) = return(x[5])

function grow_kill_fxn(n0::Int64, nS::Float64, nR::Float64, nE::Float64,
    b::Float64, d::Float64, mu::Float64, sig::Float64, del::Float64, al::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t0::Float64, tmax::Float64, t_Pass::Float64, Nmax::Int64, Cc::Int64,
    treat_ons::Array{Float64}, treat_offs::Array{Float64}, 
    Nswitch::Int64; save_at::Float64 = 0.5, treat::Bool = false, 
    n_Pass::Int64 = 2, trans_switch::Float64 = 1e+09)

    # trans_switch can be turned down to also enable the phenotypic transitions
    # to be approximated by a deterministic ODE when (nX * theta) < trans_switch
    # where nX is the number of cells in phenotype X, and theta is the 
    # transition probability out of the phenotype X compartment. Can 
    # significantly speed up the simulation for optimising/troubleshooting. 

    0 <= mu <= 1.0 || error("mu must be between 0 and 1.")
    0 <= sig <= 1.0 || error("sig must be between 0 and 1.")
    0 <= del <= 1.0 || error("del must be between 0 and 1.")
    0 <= al <= 1.0 || error("al must be between 0 and 1.")
    0 <= (al+sig) <= 1.0 || error("al and sig musn't sum to > 1.0")

    # S rates
    bS = b
    dS = d 
    # E rates
    bE = b
    dE = d

    # trade-off realised in both rates.
    # b = b * (1-δ)
    # d = d * (1-δ)
    bR = b * (1 - del)
    dR = d * (1 - del)


    # Parameter Index Guide (with initial values):

    # bS_o = 0.0,  # 1
    # dS_o = 0.0,  # 2
    # bS_j = bS,   # 3
    # dS_j = dS,   # 4
    # bR_o = 0.0,  # 5
    # dR_o = 0.0,  # 6
    # bR_j = bR,   # 7
    # dR_j = dR,   # 8
    # bE_o = 0.0,  # 9
    # dE_o = 0.0,  # 10
    # bE_j = bE,   # 11
    # dE_j = dE,   # 12
    # Dc = Dc,     # 13
    # kp = 0.0,    # 14
    # psi = psi,   # 15        
    # mu_o = 0.0,  # 16
    # sig_o = 0.0, # 17
    # al_o = 0.0,  # 18
    # mu_j = mu,   # 19
    # sig_j = sig, # 20
    # al_j = al    # 21

    p = [0.0, 0.0, bS, dS, 
         0.0, 0.0, bR, dR, 
         0.0, 0.0, bE, dE, 
         Dc, 0.0, psi,              
         0.0, 0.0, 0.0, 
         mu, sig, al] 


    u0 = [0.0, nS, nR, nE, 1]      # u[5] = Passage number
    tspan = (t0, tmax)
    
    # ODE function
    function ode_fxn!(du, u, p, t)
        
        bS_o,dS_o,bS_j,dS_j,bR_o,dR_o,bR_j,dR_j,bE_o,dE_o,bE_j,dE_j,Dc,kp,psi,mu_o,sig_o,al_o,mu_j,sig_j,al_j = p

        # Modify death rates according to the current drug concentration and 
        # total population size (re logistic growth). 

        dS_o_mod = ((p[2] + ((p[2]/dS)*p[13]*u[1])) * (1 - ((u[2]+u[3]+u[4])/Cc)))
        dR_o_mod = ((p[6] + ((p[6]/dR)*p[13]*u[1]*(1 - p[15]))) * (1 - ((u[2]+u[3]+u[4])/Cc)))
        dE_o_mod = ((p[10] + ((p[10]/dE)*p[13]*u[1]*(1 - p[15]))) * (1 - ((u[2]+u[3]+u[4])/Cc)))

        # Modify birth rates according to total population size (re logistic 
        # growth).

        bS_o_mod = (p[1]*(1 - ((u[2]+u[3]+u[4])/Cc)))
        bR_o_mod = (p[5]*(1 - ((u[2]+u[3]+u[4])/Cc)))
        bE_o_mod = (p[9]*(1 - ((u[2]+u[3]+u[4])/Cc)))

        du[1] = kp
        
        # dnS/dt
        du[2] = (bS_o_mod - dS_o_mod)*u[2] - 
                 mu_o*u[2]*bS_o_mod + sig_o*u[3]*bR_o_mod
        # dnR/dt
        du[3] = (bR_o_mod - dR_o_mod)*u[3]+ 
                 mu_o*u[2]*bS_o_mod - sig_o*u[3]*bR_o_mod - 
                 al_o*u[1]*u[3]*bR_o_mod
        # dnE/dt
        du[4] = (bE_o_mod - dE_o_mod)*u[4] + 
                 al_o*u[1]*u[3]*bR_o_mod

    end

    prob = ODEProblem(ode_fxn!, u0, tspan, p)

    # Jump Problem Birth-Death Jumps
    ################################

    # Sensitive:
    function S_birth!(integrator)
        integrator.u[2] = integrator.u[2] + 1
        nothing
    end
    function S_death!(integrator)
        integrator.u[2] = integrator.u[2] - 1
        nothing
    end
    Sb_rate(u, p, t) = (u[2] * p[3]) * (1 - ((u[2]+u[3]+u[4])/Cc))
    Sd_rate(u, p, t) = ((u[2] * p[4]) + (((p[4]/dS)*p[13]*u[1])*u[2])) * (1 - ((u[2]+u[3]+u[4])/Cc))
    Sb_jump = VariableRateJump(Sb_rate, S_birth!)
    Sd_jump = VariableRateJump(Sd_rate, S_death!)

    # Resistant:
    function R_birth!(integrator)
        integrator.u[3] = integrator.u[3] + 1
        nothing
    end
    function R_death!(integrator)
        integrator.u[3] = integrator.u[3] - 1
        nothing
    end
    Rb_rate(u, p, t) = (u[3] * p[7]) * (1 - ((u[2]+u[3]+u[4])/Cc))
    Rd_rate(u, p, t) = ((u[3] * p[8]) + (((p[8]/dR)*p[13]*u[1]*(1-p[15]))*u[3])) * (1 - ((u[2]+u[3]+u[4])/Cc))
    Rb_jump = VariableRateJump(Rb_rate, R_birth!)
    Rd_jump = VariableRateJump(Rd_rate, R_death!)

    # Escape:
    function E_birth!(integrator)
        integrator.u[4] = integrator.u[4] + 1
        nothing
    end
    function E_death!(integrator)
        integrator.u[4] = integrator.u[4] - 1
        nothing
    end
    Eb_rate(u, p, t) = (u[4] * p[11]) * (1 - ((u[2]+u[3]+u[4])/Cc))
    Ed_rate(u, p, t) = ((u[4] * p[12]) + (((p[12]/dE)*p[13]*u[1]*(1-p[15]))*u[4])) * (1 - ((u[2]+u[3]+u[4])/Cc))
    Eb_jump = VariableRateJump(Eb_rate, E_birth!)
    Ed_jump = VariableRateJump(Ed_rate, E_death!)
    
    # Jump Problem Phenotype Switching Jumps
    ########################################

    # Sensitive -> Resistant 
    function SR_switch!(integrator)
        integrator.u[2] = integrator.u[2] - 1
        integrator.u[3] = integrator.u[3] + 1
        nothing
    end

    SR_switch_rate(u, p, t) = (u[2] * p[19] * (p[1]+p[3])) * (1 - ((u[2]+u[3]+u[4])/Cc))
    SR_switch_jump = VariableRateJump(SR_switch_rate, SR_switch!)

    # Resistant -> Sensitive
    function RS_switch!(integrator)
        integrator.u[3] = integrator.u[3] - 1
        integrator.u[2] = integrator.u[2] + 1
        nothing
    end

    RS_switch_rate(u, p, t) = (u[3] * p[20] * (p[5]+p[7])) * (1 - ((u[2]+u[3]+u[4])/Cc))
    RS_switch_jump = VariableRateJump(RS_switch_rate, RS_switch!)

    # Resistant -> Escape
    function RE_switch!(integrator)
        integrator.u[3] = integrator.u[3] - 1
        integrator.u[4] = integrator.u[4] + 1
        nothing
    end

    RE_switch_rate(u, p, t) = (u[1] * u[3] * p[21] * (p[5]+p[7])) * (1 - ((u[2]+u[3]+u[4])/Cc))
    RE_switch_jump = VariableRateJump(RE_switch_rate, RE_switch!)
    

    ###########
    # Callbacks
    ###########

    # Sensitive
    ############

    # Switch from SJM -> ODE/SDE
    S_cond_switch1(u, t, integrator) = integrator.u[2] >= Nswitch

    # Affect for switching - turn the ODE integrator rates on and SJM rates off. 
    function S_switch_1!(integrator)
        if integrator.p[1] == 0.0
            integrator.p[1] = deepcopy(integrator.p[3])
            integrator.p[3] = 0.0
        end
        if integrator.p[2] == 0.0
            integrator.p[2] = deepcopy(integrator.p[4])
            integrator.p[4] = 0.0
        end
        nothing
    end

    S_cb_switch1 = DiscreteCallback(S_cond_switch1, S_switch_1!, 
                                    save_positions = (false, true))

    # Switch from ODE/SDE -> SJM
    S_cond_switch2(u, t, integrator) = integrator.u[2] < Nswitch

    # Affect for switching - turning the SJM integrator rates on and ODE rates off. 
    function S_switch_2!(integrator)
        # Round the popoulation size to discrete for the sjm.
        integrator.u[2] = round(integrator.u[2])
        if integrator.p[3] == 0.0
            integrator.p[3] = deepcopy(integrator.p[1])
            integrator.p[1] = 0.0
        end
        if integrator.p[4] == 0.0
            integrator.p[4] = deepcopy(integrator.p[2])
            integrator.p[2] = 0.0
        end
        nothing
    end

    S_cb_switch2 = DiscreteCallback(S_cond_switch2, S_switch_2!,
                                    save_positions = (false, true))

    # Resistant
    ###########

    # Switch from SJM -> ODE/SDE
    R_cond_switch1(u, t, integrator) = integrator.u[3] >= Nswitch

    # Affect for switching - turn the ODE integrator rates on and SJM rates off. 
    function R_switch_1!(integrator)
        if integrator.p[5] == 0.0
            integrator.p[5] = deepcopy(integrator.p[7])
            integrator.p[7] = 0.0
        end
        if integrator.p[6] == 0.0
            integrator.p[6] = deepcopy(integrator.p[8])
            integrator.p[8] = 0.0
        end
        nothing
    end

    R_cb_switch1 = DiscreteCallback(R_cond_switch1, R_switch_1!, 
                                    save_positions = (false, true))

    # Switch from ODE/SDE -> SJM
    R_cond_switch2(u, t, integrator) = integrator.u[3] < Nswitch

    # Affect for switching - turning the SJM integrator rates on and ODE rates off. 
    function R_switch_2!(integrator)
        # Round the popoulation size to discrete for the sjm.
        integrator.u[3] = round(integrator.u[3])
        if integrator.p[7] == 0.0
            integrator.p[7] = deepcopy(integrator.p[5])
            integrator.p[5] = 0.0
        end
        if integrator.p[8] == 0.0
            integrator.p[8] = deepcopy(integrator.p[6])
            integrator.p[6] = 0.0
        end
        nothing
    end

    R_cb_switch2 = DiscreteCallback(R_cond_switch2, R_switch_2!,
                                    save_positions = (false, true))


    # Escape
    ###########

    # Switch from SJM -> ODE/SDE
    E_cond_switch1(u, t, integrator) = integrator.u[4] >= Nswitch

    # Affect for switching - turn the ODE integrator rates on and SJM rates off. 
    function E_switch_1!(integrator)
        if integrator.p[9] == 0.0
            integrator.p[9] = deepcopy(integrator.p[11])
            integrator.p[11] = 0.0
        end
        if integrator.p[10] == 0.0
            integrator.p[10] = deepcopy(integrator.p[12])
            integrator.p[12] = 0.0
        end
        nothing
    end

    E_cb_switch1 = DiscreteCallback(E_cond_switch1, E_switch_1!, 
                                    save_positions = (false, true))

    # Switch from ODE/SDE -> SJM
    E_cond_switch2(u, t, integrator) = integrator.u[4] < Nswitch

    # Affect for switching - turning the SJM integrator rates on and ODE rates off. 
    function E_switch_2!(integrator)
        # Round the popoulation size to discrete for the sjm.
        integrator.u[4] = round(integrator.u[4])
        if integrator.p[11] == 0.0
            integrator.p[11] = deepcopy(integrator.p[9])
            integrator.p[9] = 0.0
        end
        if integrator.p[12] == 0.0
            integrator.p[12] = deepcopy(integrator.p[10])
            integrator.p[10] = 0.0
        end
        nothing
    end

    E_cb_switch2 = DiscreteCallback(E_cond_switch2, E_switch_2!,
                                    save_positions = (false, true))


    # Sensitive -> Resistant
    ########################

    # Switch from SJM -> ODE/SDE
    SR_cond_switch1(u, t, integrator) = (integrator.u[2]*mu) >= trans_switch

    # Affect for switching - turn the ODE integrator rates on and SJM rates off. 
    function SR_switch_1!(integrator)
        if integrator.p[16] == 0.0
            integrator.p[16] = deepcopy(integrator.p[19])
            integrator.p[19] = 0.0
        end
        nothing
    end

    SR_cb_switch1 = DiscreteCallback(SR_cond_switch1, SR_switch_1!, 
                                     save_positions = (false, true))

    # Switch from ODE/SDE -> SJM
    SR_cond_switch2(u, t, integrator) = (integrator.u[2]*mu) < trans_switch

    # Affect for switching - turning the SJM integrator rates on and ODE rates off. 
    function SR_switch_2!(integrator)
        if integrator.p[19] == 0.0
            integrator.p[19] = deepcopy(integrator.p[16])
            integrator.p[16] = 0.0
        end
        nothing
    end

    SR_cb_switch2 = DiscreteCallback(SR_cond_switch2, SR_switch_2!,
                                     save_positions = (false, true))

    # Resistant -> Sensitive 
    ########################

    # Switch from SJM -> ODE/SDE
    RS_cond_switch1(u, t, integrator) = (integrator.u[3]*sig) >= trans_switch

    # Affect for switching - turn the ODE integrator rates on and SJM rates off. 
    function RS_switch_1!(integrator)
        if integrator.p[17] == 0.0
            integrator.p[17] = deepcopy(integrator.p[20])
            integrator.p[20] = 0.0
        end
        nothing
    end

    RS_cb_switch1 = DiscreteCallback(RS_cond_switch1, RS_switch_1!, 
                                     save_positions = (false, true))

    # Switch from ODE/SDE -> SJM
    RS_cond_switch2(u, t, integrator) = (integrator.u[3]*sig) < trans_switch

    # Affect for switching - turning the SJM integrator rates on and ODE rates off. 
    function RS_switch_2!(integrator)
        if integrator.p[20] == 0.0
            integrator.p[20] = deepcopy(integrator.p[17])
            integrator.p[17] = 0.0
        end
        nothing
    end

    RS_cb_switch2 = DiscreteCallback(RS_cond_switch2, RS_switch_2!,
                                     save_positions = (false, true))


    # Resistant -> Escape 
    #####################

    # Switch from SJM -> ODE/SDE
    RE_cond_switch1(u, t, integrator) = (integrator.u[3]*al) >= trans_switch

    # Affect for switching - turn the ODE integrator rates on and SJM rates off. 
    function RE_switch_1!(integrator)
        if integrator.p[18] == 0.0
            integrator.p[18] = deepcopy(integrator.p[21])
            integrator.p[21] = 0.0
        end
        nothing
    end

    RE_cb_switch1 = DiscreteCallback(RE_cond_switch1, RE_switch_1!, 
                                     save_positions = (false, true))

    # Switch from ODE/SDE -> SJM
    RE_cond_switch2(u, t, integrator) = (integrator.u[3]*al) < trans_switch

    # Affect for switching - turning the SJM integrator rates on and ODE rates off. 
    function RE_switch_2!(integrator)
        if integrator.p[21] == 0.0
            integrator.p[21] = deepcopy(integrator.p[18])
            integrator.p[18] = 0.0
        end
        nothing
    end

    RE_cb_switch2 = DiscreteCallback(RE_cond_switch2, RE_switch_2!,
                                     save_positions = (false, true))


    # Treatment
    ###########

    # On
    ####

    condition1(u, t, integrator) = t ∈ treat_ons
    function affect1!(integrator)
        integrator.p[14] = k 
    end

    # Off
    #####

    condition2(u, t, integrator) = t ∈ treat_offs
    function affect2!(integrator)
        integrator.p[14] = -k
    end

    # Don't let above 1.0
    condition3a(u, t, integrator) = (u[1] - 1.0)
    function affect3a!(integrator)
        integrator.u[1] = 1.0  
        integrator.p[14] = 0.0 
    end

    # Don't let below 0.0
    condition3b(u, t, integrator) = (u[1])
    function affect3b!(integrator)
        integrator.u[1] = 0.0  
        integrator.p[14] = 0.0 
    end

    cb1 = DiscreteCallback(condition1,affect1!, 
                            save_positions=(false, true))
    cb2 = DiscreteCallback(condition2,affect2!, 
                            save_positions=(false, true))
    cb3a = ContinuousCallback(condition3a,affect3a!, 
                            save_positions=(false, true))
    cb3b = ContinuousCallback(condition3b,affect3b!, 
                            save_positions=(false, true))

    # Nmax, Extinction & Passage
    ############################

    # Nmax
    ######
    # If carrying capacity reached, either passage (if current < n_Pass) or
    # terminate. 

    condition4(u, t, integrator) = ((integrator.u[2] + integrator.u[3] + integrator.u[4]) >= Nmax)

    function affect4!(integrator)
        # If currently already on n_Pass, terminate integrator
        if integrator.u[5] >= n_Pass
            terminate!(integrator)
        else
            # Sample cells without replacement into next Passage
            curr_nR = Int64(round(integrator.u[3]))
            curr_nS = Int64(round(integrator.u[2]))
            curr_nE = Int64(round(integrator.u[4]))

            exp_nRs = repeat(["R"], curr_nR)
            exp_nSs = repeat(["S"], curr_nS)
            exp_nEs = repeat(["E"], curr_nE)

            exp_ns = vcat(exp_nRs, exp_nSs, exp_nEs)
            samp_ns = sample(exp_ns, n0, replace=false)

            samp_nR = sum(samp_ns .== "R")
            samp_nS = sum(samp_ns .== "S")
            samp_nE = sum(samp_ns .== "E")

            integrator.u[3] = Float64(samp_nR)
            integrator.u[2] = Float64(samp_nS)
            integrator.u[4] = Float64(samp_nE)

            # Update the passage tracker 
            integrator.u[5] += 1
        end
    end

    cb4 = DiscreteCallback(condition4, affect4!,
                        save_positions=(false,true))

    # Extinction 
    ############
    # End the simulation if extinction occurs (both n sum to < 1.0). 

    condition5(u, t, integrator) = ((integrator.u[2] + integrator.u[3] + integrator.u[4]) < 1.0)

    cb5 = DiscreteCallback(condition5, terminate!,
                        save_positions=(false,false))
    
    # Passage time
    ##############
    # Passage cells if t_Pass has been reached and sufficient cells. Otherwise, 
    # continue growing until tmax, Nmax or extinction. 

    condition6(u, t, integrator) = t == t_Pass

    function affect6!(integrator)

        # Only passage if current Passage < n_Pass
        if integrator.u[5] < n_Pass
            # Sample cells without replacement into next Passage
            curr_nR = round(integrator.u[3])
            curr_nS = round(integrator.u[2])
            curr_nE = round(integrator.u[4])
            # Only sample if enough cells (>= n0). Otherwise, keep growing. 
            if (curr_nR+curr_nS+curr_nE) >= n0
                # Sample cells without replacement into next Passage
                curr_nR = Int64(round(integrator.u[3]))
                curr_nS = Int64(round(integrator.u[2]))
                curr_nE = Int64(round(integrator.u[4]))
    
                exp_nRs = repeat(["R"], curr_nR)
                exp_nSs = repeat(["S"], curr_nS)
                exp_nEs = repeat(["E"], curr_nE)

                exp_ns = vcat(exp_nRs, exp_nSs, exp_nEs)
                samp_ns = sample(exp_ns, n0, replace=false)

                samp_nR = sum(samp_ns .== "R")
                samp_nS = sum(samp_ns .== "S")
                samp_nE = sum(samp_ns .== "E")

                integrator.u[3] = Float64(samp_nR)
                integrator.u[2] = Float64(samp_nS)
                integrator.u[4] = Float64(samp_nE)
    
                # Update the passage tracker 
                integrator.u[5] += 1
            end
        end
    end

    cb6 = DiscreteCallback(condition6, affect6!, 
                           save_positions = (false, true))
                
        

    # Population Limits
    ###################

    # Don't let below 0.0
    condition7a(u, t, integrator) = u[2] < 0
    function affect7a!(integrator)
        integrator.u[2] = 0.0  
    end
    cb7a = DiscreteCallback(condition7a,affect7a!, 
                              save_positions=(false, true))

    condition7b(u, t, integrator) = u[3] < 0
    function affect7b!(integrator)
        integrator.u[3] = 0.0  
    end
    cb7b = DiscreteCallback(condition7b,affect7b!, 
                            save_positions=(false, true))

    condition7c(u, t, integrator) = u[4] < 0
    function affect7c!(integrator)
        integrator.u[4] = 0.0  
    end
    cb7c = DiscreteCallback(condition7c,affect7c!, 
                                save_positions=(false, true))
                        

    # Collect Callbacks 
    ###################

    # Callbacks into a single set - include treatment if treat == true
    if treat == true
        cbs = CallbackSet(S_cb_switch1, S_cb_switch2, 
                          R_cb_switch1, R_cb_switch2,
                          E_cb_switch1, E_cb_switch2, 
                          SR_cb_switch1, SR_cb_switch2, 
                          RS_cb_switch1, RS_cb_switch2,
                          RE_cb_switch1, RE_cb_switch2, 
                          cb1, cb2, cb3a, cb3b, cb4, cb5, cb6,
                          cb7a, cb7b, cb7c)
    else
        cbs = CallbackSet(S_cb_switch1, S_cb_switch2, 
                          R_cb_switch1, R_cb_switch2,
                          E_cb_switch1, E_cb_switch2, 
                          SR_cb_switch1, SR_cb_switch2, 
                          RS_cb_switch1, RS_cb_switch2,
                          RE_cb_switch1, RE_cb_switch2, 
                          cb4, cb5, cb6,
                          cb7a, cb7b, cb7c)
    end

    #######
    # Solve
    #######

    # Convert into a Jump problem
    sjm_prob = JumpProblem(prob, Direct(), 
                           Sb_jump, Sd_jump, 
                           Rb_jump, Rd_jump,
                           Eb_jump, Ed_jump,
                           SR_switch_jump, 
                           RS_switch_jump,
                           RE_switch_jump)

    sol = @suppress solve(sjm_prob, Tsit5(), callback = cbs, 
    tstops = collect(0:save_at:tmax))
    
    return sol

end


# Run the growth-kill function multiple times for each replicate, and only 
# return the population sizes and respective times used for parameter 
# estimation. 

function exp_pulse_treat_fxn(n0::Int64, b::Float64, d::Float64,
    rho::Float64, mu::Float64, sig::Float64, del::Float64, al::Float64,
    Dc::Float64, k::Float64, psi::Float64,
    t0::Float64, t_exp::Float64, tmax::Float64, t_Pass::Float64, 
    Nmax::Int64, Cc::Int64,
    treat_ons::Array{Float64}, treat_offs::Array{Float64}, 
    t_keep::Array{Float64}, Nswitch::Int64; 
    save_at::Float64 = 0.5, n_Pass::Int64 = 2)

    # If sig and al sum to > 1.0, throw some dummy data that will produce
    # a large distance in the inference: 

    if (al+sig) > 1.0

        fin_t_outs = repeat([-1.0], 4*(length(t_keep)+n_Pass))
        fin_u_outs = repeat([-1], 4*(length(t_keep)+n_Pass))

    else

        # Calculate nS and nR using n0 and rho. 
        nR = Float64(floor(rho * n0))
        nS = n0 - nR
        # nE = always 0.0 at t = 0.0. 
        nE = 0.0

        # Expand the cells in an untrated environment: 
        exp_cells = grow_kill_fxn(n0, nS, nR, nE, b, d, mu, sig, del, al, 
                                Dc, k, psi, 0.0, 
                                t_exp, -1.0, Nmax, Cc,
                                treat_ons, treat_offs, Nswitch, 
                                save_at=save_at, treat=false, n_Pass=1)

        # Extract number of sensitive and resistant cells: 
        exp_nS = Int64(round(pop_fun(last(exp_cells.u))[1]))
        exp_nR = Int64(round(pop_fun(last(exp_cells.u))[2]))
        exp_nE = Int64(round(pop_fun(last(exp_cells.u))[3]))

        # Check enough cells for sampling. If not, throw some dummy data that we 
        # can use to throw an impossibly high distance in the ABC step. 
        if (exp_nS + exp_nR + exp_nE) < n0*4

            # Populate output arrays with negative values (impossible in the 
            # true sim). 
            fin_t_outs = repeat([-1], 4*(length(t_keep)+n_Pass))
            fin_u_outs = repeat([-1], 4*(length(t_keep)+n_Pass))
            # Return the vectors of t and u
            return Dict("t"=>fin_t_outs, "u"=>fin_u_outs)

        else
        
            # Create matrix to hold sampled phenotypes per replicate
            rep_phenos = Array{Float64}(undef, 3, 4)

            exp_nSs = repeat(["S"], exp_nS)
            exp_nRs = repeat(["R"], exp_nR)
            exp_nEs = repeat(["E"], exp_nE)

            exp_ns = vcat(exp_nRs, exp_nSs, exp_nEs)

            # Go through and sample nRs, removing from exp_nR each time. 
            for i in 1:4
                # Sample without replacement 
                samp_ids = sample(1:length(exp_ns), n0, replace=false)
                samp_ns = exp_ns[samp_ids]

                samp_nS = sum(samp_ns .== "S")
                samp_nR = sum(samp_ns .== "R")
                samp_nE = sum(samp_ns .== "E")

                # Add to replicate starting phenotype proportions
                rep_phenos[1,i] = Float64(samp_nS)
                rep_phenos[2,i] = Float64(samp_nR)
                rep_phenos[3,i] = Float64(samp_nE)

                # Remove the sample ids from the exp_n vector
                deleteat!(exp_ns, sort(samp_ids))
            end

            # Initialise start time 
            t0 = 0.0

            # Vectors to keep final outputs: 
            fin_t_outs = Float64[]
            fin_u_outs = Int64[]

            # Repeat simulation for each replicate:
            for i in 1:4

                # Run to derive solution 
                sol = grow_kill_fxn(n0, 
                                    rep_phenos[1,i], 
                                    rep_phenos[2,i], 
                                    rep_phenos[3,i],
                                    b, d, mu, sig, del, al, Dc, k, psi, t0, 
                                    tmax, t_Pass, Nmax, Cc,
                                    treat_ons, treat_offs,
                                    Nswitch, save_at=save_at, treat=true, 
                                    n_Pass=n_Pass)

                # Vectors to save output t and population sizes (u) - size of 
                # t_keep vector, plus number of passages. 
                t_outs = Vector{Float64}(undef, length(t_keep)+n_Pass)
                u_outs = Vector{Int64}(undef, length(t_keep)+n_Pass)


                # Save the t_keep times. 
                for j in 1:length(t_keep)
                    # If the current time wasn't in the solution, find closest 
                    # value to keep. 
                    if !(t_keep[j] in sol.t)
                        # Find the closest time reached for this value in t_keep
                        t_closest_pos = findmin(abs.(sol.t .- t_keep[j]))[2]
                        t_closest = sol.t[t_closest_pos]
                        t_outs[j+2] = t_closest
                        # Keep the pop size reached at this time
                        u_outs[j+2] = Int64(round(sum.(pop_fun.(sol.u))[t_closest_pos]))
                    # Otherwise keep the closest solution 
                    else
                        t_realised_pos = findlast(sol.t .== t_keep[j][1])
                        t_realised = sol.t[t_realised_pos]
                        t_outs[j+2] = t_realised
                        u_outs[j+2] = Int64(round(sum.(pop_fun.(sol.u))[t_realised_pos]))
                    end
                end

                # Extract the Passage times using the given integrator variable. 
                # If Passage not reached, use the most recent time. 
                for j in 1:n_Pass
                    if sum(pass_fun.(sol.u) .== j) > 0
                        t_realised_pos = findlast(pass_fun.(sol.u) .== j)
                        t_outs[length(t_keep)+j] = sol.t[t_realised_pos]
                        u_outs[length(t_keep)+j] = Int64(round(sum.(pop_fun.(sol.u))[t_realised_pos]))
                    else
                        # If current passage not realised, save previous pop size. 
                        t_outs[length(t_keep)+j] = t_outs[length(t_keep)+j-1]
                        u_outs[length(t_keep)+j] = u_outs[length(t_keep)+j-1]

                    end
                end
                append!(fin_t_outs, t_outs)
                append!(fin_u_outs, u_outs)
            end

            # Round t
            fin_t_outs = Float64.(round.(fin_t_outs, digits = 0))

        end

    end

    # Return the vectors of t and u
    return Dict("t"=>fin_t_outs, "u"=>fin_u_outs)
    
end


# Measurement Noise
###################

function model_meas_noise(u_outs::Array{Int64}, phi::Float64)

    # If any negatives values exist pre-noise, just return original values. 
    if sum(u_outs .< 0) > 0

        noise_u_outs = u_outs

    # Otherwise calculate distance normally.
    else

        # Sample the numbers + noise (where sd = phi * obsv_N)

        noise_u_outs = rand.(Normal.(u_outs, u_outs*phi))

        # Replace any negatives with 0s. 

        noise_u_outs[noise_u_outs .< 0] .= 0

        # Round
        noise_u_outs = Int64.(round.(noise_u_outs, digits = -2))

    end

    return noise_u_outs

end


################################################################################

########################################################
# Set the fixed parameters for the simualted experiment:
########################################################

# Number of uniquely barcoded cells when the experiment begins:
n0 = Int64(1e+04);

# Estimates of the average, sensitive population's birth and death rates:
b = 0.893; d = 0.200;

# The measurement error of the population size readings (where
# observed N ~ Normal(mean = true N, sd = true N * phi)):
phi = 0.10;

# Time when the experiment begins and maximum time of the experiment (in days):
t0 = 0.0; tmax = 60.0;

# Maximum number of cells before a replicate is harvested, and an estimate of
# the carrying capacity:
Nmax = Int64(80*1e+04); Cc = Int64(100*1e+04);

# The treatment windows when treatment begins (treat_ons) and ends (treat_offs)
# (in days):
treat_ons = collect(4.0:8.0:tmax);
treat_offs = collect(8.0:8.0:tmax);

# The time to Passage cells (in days) and number of Passages (currently 
# max n_Pass = 2):
t_Pass = 30.0; n_Pass = 2

# The expansion time of the barcoded cells before splitting into experimental
# replicates (in days):
t_exp = 6.0;

# The observation times to record the population size changes (in days):
t_keep = [5.0, 10.0];

# The time-window when calculating the drug-concentration change, and a Boolean
# whether to treat the cells (leave these unless troubleshooting):
dt_save_at = 1e-03; treat = true; 

# The population size per-phenotype to switch from a stochastic jump process
# (when <= Nswitch) to a deterministic ODE approximation (> Nswitch).
Nswitch = 10;

################################################################################



"""
Simulate model for given parameters. 
"""

function model(pars)

    # Choose which model to run based on the parameters provided: 

    # C - Escape Transitions
    ########################

    if "al" in keys(pars)

        if pars["sig"] == -1
            sig = 0.0
        else
            sig = 10.0.^((-pars["sig"]))
        end

        if pars["al"] == -1
            al = 0.0
        else
            al = 10.0.^((-pars["al"]))
        end

        # Extract parameters and transform off the log scale. 
        rho = 10.0.^((-pars["rho"]))
        mu = 10.0.^((-pars["mu"]))
        Dc = convert(Float64, pars["Dc"])
        del = 10.0.^((-pars["delt"]))
        k = convert(Float64, pars["k"])
        psi = convert(Float64, pars["psi"])

        # Run the simulation
        sim_out = exp_pulse_treat_fxn(n0, b, d, rho, mu, sig, del, al, 
                                      Dc, k, psi,
                                      t0, t_exp, tmax, t_Pass, 
                                      Nmax, Cc,
                                      treat_ons, treat_offs, t_keep, Nswitch)

    # B - Bidirectional Transitions
    ###############################

    elseif "sig" in keys(pars)

        if pars["sig"] == -1
            sig = 0.0
        else
            sig = 10.0.^((-pars["sig"]))
        end

        # Extract parameters and transform off the log scale. 
        rho = 10.0.^((-pars["rho"]))
        mu = 10.0.^((-pars["mu"]))
        Dc = convert(Float64, pars["Dc"])
        del = 10.0.^((-pars["delt"]))
        k = convert(Float64, pars["k"])
        psi = convert(Float64, pars["psi"])

        al = 0.0

        # Run the simulation
        sim_out = exp_pulse_treat_fxn(n0, b, d, rho, mu, sig, del, al, 
                                      Dc, k, psi,
                                      t0, t_exp, tmax, t_Pass, 
                                      Nmax, Cc,
                                      treat_ons, treat_offs, t_keep, Nswitch)


    # A - Unidirectional Transitions
    ################################

    else

        sig = 0.0
        al = 0.0

        # Extract parameters and transform off the log scale. 
        rho = 10.0.^((-pars["rho"]))
        mu = 10.0.^((-pars["mu"]))
        Dc = convert(Float64, pars["Dc"])
        del = 10.0.^((-pars["delt"]))
        k = convert(Float64, pars["k"])
        psi = convert(Float64, pars["psi"])

        # Run the simulation
        sim_out = exp_pulse_treat_fxn(n0, b, d, rho, mu, sig, del, al, 
                                      Dc, k, psi,
                                      t0, t_exp, tmax, t_Pass, 
                                      Nmax, Cc,
                                      treat_ons, treat_offs, t_keep, Nswitch)

    end

    # Add measurement noise
    sim_out["u"] = model_meas_noise(sim_out["u"], phi)

    # Return output as dictionary 
    return Dict("t"=>sim_out["t"], "u"=>sim_out["u"])
end


"""
2D Euclidean Distance between model simulations or observed data `y` and `y0`.
Uses vector of times (in hrs) and log10(Nt + 1) to calcaulte distances, 
using the replicate permutation that generates the lowest distance.
"""

function distance(y, y0)

    t, t0 = y["t"], y0["t"]
    u, u0 = y["u"], y0["u"]

    # Split the vectors up into their replicates' values:

    function split_vector(vec, n)
        len = length(vec)
        # Ensure that the vector can be divided into n equal parts
        if len % n != 0
            error("The length of the vector is not divisible by n")
        end

        part_size = len ÷ n  # Integer division
        return [vec[i:i+part_size-1] for i in 1:part_size:len]
    end

    split_t0 = split_vector(t0, 4)
    split_u0 = split_vector(u0, 4)
    
    split_t = split_vector(t, 4)
    split_u = split_vector(u, 4)
    
    # Calculate the average Euclidean distance given a single 
    # replicate's t and u values: 

    function euc_dist(t0, u0, t, u)

        # Check all vectors the same length.
        length(t0) == length(u0) == length(t) == length(u) || error("Dimension mismatch.")

        # If negative values in any of the vectors, throw a v.high distance. 
        if sum(vcat(t .< 0, t0 .< 0)) > 0 || sum(vcat(u .< 0, u0 .< 0)) > 0

            D = 10^10

        # Otherwise calculate distance normally.
        else

            # Put pop sizes on log scale
            u, u0 = log10.(u.+1), log10.(u0.+1)

            # Normalise time and pop vectors by tmax and Nmax
            u, u0 = u/(log10(Nmax+1)), u0/(log10(Nmax+1))

            t, t0 = t/(tmax*1.0), t0/(tmax*1.0)

            dist = sqrt.(((t .- t0).^2) + ((u .- u0).^2))

            D = mean(dist)

        end

        return D

    end

    # Now calculate the mean distance for every iteration of the simulation 
    # replicates:

    sequence = [1, 2, 3, 4]
    all_permutations = []

    for p in permutations(sequence, 4)
        push!(all_permutations, p)
    end

    # Now, given the observed replicates time (t) and population (u) vectors, 
    # calculate the average distance across all 4 replicates given all the 
    # possible permutations of the simualted replicates:

    perm_dist = zeros(length(all_permutations))

    for i in eachindex(all_permutations)

        temp_rep_dist = 0.0

        for j in 1:4
            
            rep_dist = euc_dist(split_t0[j], 
                                split_u0[j],
                                split_t[all_permutations[i][j]],
                                split_u[all_permutations[i][j]])

            temp_rep_dist += rep_dist

        end

        # Need to nromalise by number of replicates:
        perm_dist[i] = temp_rep_dist/4

    end

    # Return the minimum distance 

    min_D = minimum(perm_dist)

    return(min_D)

end

end  # module


