module BlkTX
using UnPack, DSP, Random, Interpolations
include("../util/Util_JLSD.jl"); using .Util_JLSD

export dac_drv_top!

function dac_drv_top!(drv, Si)
    @unpack osr, dt, blk_size, blk_size_osr = drv.param
    @unpack ir, fir_norm, swing, Sfir_mem, Vo_mem = drv

    u_filt!(drv.Sfir_conv, Si, fir_norm, Si_mem = Sfir_mem)
    Vi_fir = kron(drv.Sfir, ones(osr))

    
    if drv.jitter_en
        drv_add_jitter!(drv, Vi_fir)
    end

    u_conv!(drv.Vo_conv, Vi_fir, ir, Vi_mem=Vo_mem, gain=dt*swing/2)

end

function drv_add_jitter!(drv, Vi_fir)
    @unpack osr, tui, blk_size = drv.param
    @unpack dcd, rj_s, last_rj_osr, sj_amp_ui, sj_freq, last_sj_phi = drv
    @unpack prev_nui, V_prev_nui = drv

    rj_osr = rj_s/tui*osr
    sj_amp_osr = sj_amp_ui*osr
    sj_freq_osr = sj_freq*tui/osr

    tvec, drv.last_rj_osr, drv.last_sj_phi = drv_jitter_tvec(blk_size, osr, 
                                        ;dcd, rj_osr, last_rj_osr,
                                        sj_amp_osr, sj_freq_osr, last_sj_phi)
    drv.Vext .= [V_prev_nui; Vi_fir]
    drv.V_prev_nui .= Vi_fir[end-prev_nui*osr+1:end]

    Vi_fir .= linear_interpolation(drv.tt_Vext, drv.Vext).(tvec)
end


function drv_jitter_tvec(len, osr; dcd, rj_osr, last_rj_osr, sj_amp_osr, sj_freq_osr, last_sj_phi)

    t_ui::Vector{Float64} = LinRange(0,len*osr,len+1)
    t_rj::Vector{Float64} = [last_rj_osr; rj_osr*randn(len)]
    phi_sj::Vector{Float64} = (last_sj_phi .+ 2*π*sj_freq_osr .* t_ui).%(2*π)
    t_sj::Vector{Float64} = sj_amp_osr*sin.(phi_sj)

    t_ui[2:2:end] .+= osr*dcd
    t_ui .+= (t_rj .+ t_sj)

    tvec_jitter = zeros(len*osr)
    for n = 1:lastindex(t_ui)-1
        tvec_jitter[(n-1)*osr+1:n*osr] .= LinRange(t_ui[n], t_ui[n+1], osr+1)[1:end-1]
    end

    return tvec_jitter, t_rj[end], phi_sj[end]
end

end