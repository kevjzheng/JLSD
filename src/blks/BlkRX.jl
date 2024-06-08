module BlkRX
using UnPack, DSP, Random, Interpolations
include("../util/Util_JLSD.jl"); using .Util_JLSD

export clkgen_pi_itp_top!
export sample_itp_top!, sample_phi_top!, slicers_top!
export cdr_top!, adpt_top!


function clkgen_pi_itp_top!(clkgen; pi_code)
    @unpack tui, osr, cur_subblk, subblk_size = clkgen.params
    @unpack nphases, rj, skews = clkgen
    @unpack pi_code_prev, pi_wrap_ui, pi_wrap_ui_Δcode = clkgen
    @unpack pi_nonlin_lut, pi_ui_cover, pi_codes_per_ui = clkgen
    
    Δpi_code = pi_code-pi_code_prev
    if abs(Δpi_code) > pi_wrap_ui_Δcode
        pi_wrap_ui -= sign(Δpi_code)*pi_ui_cover
    end

    Φ0 = osr*(pi_wrap_ui + (pi_code + pi_nonlin_lut[pi_code+1])/pi_codes_per_ui)
    Φstart = (cur_subblk-1)*subblk_size*osr
    Φnom = Φstart:osr:Φstart+(subblk_size-1)*osr
    Φskew = kron(ones(Int(subblk_size/nphases)), skews/tui*osr)
    Φrj = rj/tui*osr*randn(subblk_size)

    @. clkgen.Φo = Φ0 + Φnom + Φskew + Φrj
    
    clkgen.pi_code_prev = pi_code
    clkgen.pi_wrap_ui = pi_wrap_ui
    
end



function sample_itp_top!(splr, Vi)
    @unpack osr,dt, blk_size_osr = splr.params
    @unpack prev_nui, V_prev_nui, ir, Vo_mem = splr
    
    u_conv!(splr.Vo, splr.Vo_mem, Vi, ir, dt, blk_size_osr, Vi_mem=Vo_mem)

    splr.Vext = [V_prev_nui; splr.Vo]
    splr.V_prev_nui = splr.Vo[end-prev_nui*osr+1:end]
    splr.itp_Vext = linear_interpolation(splr.tt_Vext, splr.Vext)
end

function sample_phi_top!(splr, Φi)
    @unpack cur_subblk, subblk_size = splr.params
    @unpack itp_Vext = splr

    splr.So_subblk = itp_Vext.(Φi)
    append!(splr.So, splr.So_subblk)

end

function slicers_top!(slc, Si; ref_code)
    @unpack nphases, noise_rms, dac_min, dac_lsb = slc
    @unpack ofsts, N_per_phi = slc

    for n = eachindex(Si)
        phi_idx = (n-1)%nphases + 1
        nslc = N_per_phi[phi_idx]

        if nslc != 0
            slc.So[n] = (Si[n]  .- (dac_min.+dac_lsb*ref_code[phi_idx])
                                .+ ofsts[phi_idx] 
                                .+ noise_rms*randn(nslc)) .> 0
        end
    end
end


function cdr_top!(cdr, Sd, Se)
    @unpack Neslc_per_phi, Sd_prev = cdr
    @unpack eslc_nvec, filt_patterns, kp, ki = cdr
    @unpack pd_accum, ki_accum, pd_gain, pi_res = cdr

    pi_bnd = 2^pi_res
    Sd_val = [Sd_prev; [sum(dvec) for dvec in Sd]]

    for n = findall(eslc_nvec.!=0)
        if Sd_val[n:n+2] in filt_patterns
            vote = sign(Se[n][1].-0.5)*sign(Sd_val[n]-Sd_val[n+2])
            ki_accum += ki*vote
            pd_accum += pd_gain*(kp*vote + ki_accum)
        end
    end

    cdr.pd_accum = (pd_accum < 0) ? pi_bnd + pd_accum : (pd_accum >= pi_bnd) ? pd_accum - pi_bnd : pd_accum
    cdr.pi_code = Int(floor(cdr.pd_accum))


    cdr.Sd_prev = Sd_val[end]    
    cdr.ki_accum = ki_accum

end

function adpt_top!(adpt, Sd, Se)
    @unpack Neslc_per_phi, Sd_prev = adpt
    @unpack eslc_nvec, eslc_filt_patterns, eslc_ref_max, mu_eslc = adpt

    Sd_val = [Sd_prev; [sum(dvec) for dvec in Sd]]

    ref_accum = adpt.eslc_ref_accum

    for n = findall(eslc_nvec.!=0)
        ref_accum +=  (Sd_val[n:n+2] in eslc_filt_patterns) ? 
                        mu_eslc*sign(Se[n][1].-0.5) : 0
    end
    adpt.eslc_ref_accum =   ref_accum < 0 ? 0 : 
                            ref_accum > eslc_ref_max ? eslc_ref_max :
                            ref_accum
    adpt.eslc_ref_code = floor(adpt.eslc_ref_accum)
    adpt.eslc_ref_vec = [adpt.eslc_ref_code*ones(Int,n) for n in Neslc_per_phi]

    adpt.Sd_prev = Sd_val[end]    
end

end