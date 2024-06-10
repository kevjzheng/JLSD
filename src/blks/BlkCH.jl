module BlkCH
using UnPack, DSP, Random
include("../util/Util_JLSD.jl"); using .Util_JLSD

export ch_top!

function ch_top!(ch,Vi)
    @unpack osr, blk_size, dt, blk_size_osr = ch.param
    @unpack noise_en, noise_rms= ch
    @unpack ir_ch, ir_pad, Vch_mem, Vo_mem = ch


    u_conv!(ch.Vch_conv, Vi, ir_ch, Vi_mem=Vch_mem, gain = dt)

    if noise_en
        ch.Vch .+= noise_rms .* randn(blk_size_osr)
    end

    u_conv!(ch.Vo_conv, ch.Vch, ir_pad, Vi_mem=Vo_mem, gain = dt)

end



end