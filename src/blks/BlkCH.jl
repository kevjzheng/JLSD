module BlkCH
using UnPack, DSP, Random
include("../util/Util_JLSD.jl"); using .Util_JLSD

export ch_top!

function ch_top!(ch,Vi)
    @unpack osr, blk_size, dt, blk_size_osr, pr_sinc = ch.params
    @unpack noise_en, noise_rms, No_mem = ch
    @unpack ir_ch, ir_pad, Vch_mem, Vo_mem = ch


    u_conv!(ch.Vch, ch.Vch_mem, Vi, ir_ch, dt, blk_size_osr, Vi_mem=Vch_mem)

    if noise_en
        noise_samples = kron(noise_rms*randn(blk_size),[1;zeros(osr-1)])
        u_conv!(ch.No, ch.No_mem, noise_samples, pr_sinc, 1, blk_size_osr, Vi_mem=No_mem)
    end

    u_conv!(ch.Vo, ch.Vo_mem, ch.Vch+ch.No, ir_pad, dt, blk_size_osr, Vi_mem=Vo_mem)

end



end