module BlkTX
using UnPack, DSP, Random, Interpolations
include("../util/Util_JLSD.jl"); using .Util_JLSD

export dac_drv_top!

function dac_drv_top!(drv, Si)
    @unpack osr, dt, blk_size, blk_size_osr = drv.param
    @unpack ir, fir_norm, swing, Vfir, Sfir_mem, Vo_mem = drv

    u_filt!(drv.Sfir_conv, Si, fir_norm, Si_mem = Sfir_mem)

    kron!(drv.Vfir, drv.Sfir, ones(osr))

    if drv.jitter_en
        drv_add_jitter!(drv, drv.Vfir)
    end

    u_conv!(drv.Vo_conv, drv.Vfir, ir, Vi_mem=Vo_mem, gain=dt*swing/2)


end

function drv_add_jitter!(drv, Vfir)
    @unpack osr, tui, blk_size, blk_size_osr = drv.param
    @unpack dcd, rj_s, sj_amp_ui, sj_freq, last_sj_phi = drv
    @unpack Vext, V_prev_nui = drv
    @unpack Δtt_ext, Δtt_prev_nui, Δtt, tt_Vext, tt_uniform = drv

    rj_osr = rj_s/tui*osr
    sj_amp_osr = sj_amp_ui*osr
    sj_freq_norm = sj_freq*tui

    drv.last_sj_phi = drv_jitter_Δt!(Δtt; blk_size, osr, 
                                        dcd, rj_osr, 
                                        sj_amp_osr, sj_freq_norm, last_sj_phi)

    Vext[eachindex(V_prev_nui)] .= V_prev_nui
    Vext[lastindex(V_prev_nui)+1:end] .= Vfir
    Δtt_ext[eachindex(Δtt_prev_nui)] .= Δtt_prev_nui
    Δtt_ext[lastindex(Δtt_prev_nui)+1:end] .= Δtt

    drv_jitter_tvec!(tt_Vext, Δtt_ext, osr)

    drv_interp_jitter!(Vfir, tt_Vext, Vext, tt_uniform)
    return nothing
end


function drv_jitter_Δt!(Δtt; blk_size, osr, dcd, rj_osr, sj_amp_osr, sj_freq_norm, last_sj_phi)

    fill!(Δtt, zero(Float64))
    Δtt[1:2:end] .+= dcd/2*osr
    Δtt[2:2:end] .-= dcd/2*osr
    Δtt .+= rj_osr .* randn(blk_size) #add rj
    
    phi_sj = (last_sj_phi .+ (2π*sj_freq_norm) * (1:blk_size)) .% (2π)
    Δtt .+= sj_amp_osr .* sin.(phi_sj) #add sj


    return phi_sj[end]
end

function drv_jitter_tvec!(tt_Vext, Δtt_ext, osr)
    for n = 1:lastindex(Δtt_ext)-1
        tt_Vext[(n-1)*osr+1:n*osr] .= LinRange((n-1)*osr+Δtt_ext[n], n*osr+Δtt_ext[n+1], osr+1)[1:end-1]
    end

    return nothing
end


function drv_interp_jitter!(vo, tt_jitter, vi, tt_uniform)
    last_idx = 1
    for n = eachindex(tt_uniform)
        t = tt_uniform[n]
        for m = last_idx:lastindex(tt_jitter)-1
            if (t >= tt_jitter[m]) && (t < tt_jitter[m+1])
                k = (vi[m+1]-vi[m])/(tt_jitter[m+1]-tt_jitter[m])
                vo[n] = vi[m] + k*(t-tt_jitter[m])
                last_idx = m
                break
            end
        end
    end

    return nothing
end

end