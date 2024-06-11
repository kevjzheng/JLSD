module BlkTX
using UnPack, DSP, Random, Interpolations
include("../util/Util_JLSD.jl"); using .Util_JLSD

export dac_drv_top!

function dac_drv_top!(drv, Si)
    @unpack osr, dt, blk_size, blk_size_osr = drv.param
    @unpack ir, fir_norm, swing, Vfir, Sfir_mem, Vo_mem = drv

    u_filt!(drv.Sfir_conv, Si, fir_norm, Si_mem = Sfir_mem)
    Vfir .= kron(drv.Sfir, ones(osr))

    
    if drv.jitter_en
        drv_add_jitter!(drv, Vfir)
    end

    u_conv!(drv.Vo_conv, Vfir, ir, Vi_mem=Vo_mem, gain=dt*swing/2)

end

function drv_add_jitter!(drv, Vfir)
    @unpack osr, tui, blk_size = drv.param
    @unpack dcd, rj_s, last_rj_osr, sj_amp_ui, sj_freq, last_sj_phi = drv
    @unpack V_prev_nui, Vext, tt_Vext, tt_jitter = drv

    rj_osr = rj_s/tui*osr
    sj_amp_osr = sj_amp_ui*osr
    sj_freq_osr = sj_freq*tui/osr

    drv.last_rj_osr, drv.last_sj_phi = drv_jitter_tvec!(tt_jitter; blk_size, osr, 
                                        dcd, rj_osr, last_rj_osr,
                                        sj_amp_osr, sj_freq_osr, last_sj_phi)
    Vext[eachindex(V_prev_nui)] .= V_prev_nui
    Vext[lastindex(V_prev_nui)+1:end] .= Vfir

    Vfir .= linear_interpolation(tt_Vext, Vext).(tt_jitter)

    return nothing
end


function drv_jitter_tvec!(tt_jitter; blk_size, osr, dcd, rj_osr, last_rj_osr, sj_amp_osr, sj_freq_osr, last_sj_phi)

    tt_ui::Vector{Float64} = LinRange(0, blk_size*osr, blk_size+1)
    tt_rj::Vector{Float64} = [last_rj_osr; rj_osr*randn(blk_size)]
    phi_sj::Vector{Float64} = (last_sj_phi .+ 2*π*sj_freq_osr .* tt_ui).%(2*π)
    tt_sj::Vector{Float64} = sj_amp_osr*sin.(phi_sj)

    tt_ui[2:2:end] .+= dcd*osr
    tt_ui .+= tt_rj .+ tt_sj

    for n = 1:lastindex(tt_ui)-1
        tt_jitter[(n-1)*osr+1:n*osr] .= LinRange(tt_ui[n], tt_ui[n+1], osr+1)[1:end-1]
    end

    return tt_rj[end], phi_sj[end]
end

end