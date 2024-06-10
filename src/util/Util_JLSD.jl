module Util_JLSD
using StatsBase, DSP, Interpolations, FFTW
using MAT

export u_conv!, u_filt!
export u_gen_ir_rc, u_fr_to_imp, u_hist


function u_conv(input, ir, dt, len; Vi_mem = zeros(1), gain = 1)
    vconv = gain*dt*conv(ir, input)
    vconv[eachindex(Vi_mem)] += Vi_mem

    return vconv[1:len], vconv[len+1:end]
end

function deprecate_u_conv!(Vo, Vo_mem, input, ir, dt, len; Vi_mem = Float64[], gain = 1)
    vconv = gain*dt*conv(input, ir)
    vconv[eachindex(Vi_mem)] .+= Vi_mem
    Vo .= vconv[1:len]
    Vo_mem .= vconv[len+1:end]

    return nothing
end

function u_conv!(Vo_conv, input, ir; Vi_mem = Float64[], gain = 1)
    Vo_conv[eachindex(Vi_mem)] .= Vi_mem
    Vo_conv[lastindex(Vi_mem)+1:end] .= zero(Float64)

    Vo_conv .+= conv(gain .* input, ir)
    return nothing
end



function u_filt!(So_conv, input, fir; Si_mem)
    So_conv[eachindex(Si_mem)] .= Si_mem
    So_conv[lastindex(Si_mem)+1:end] .= zero(Float64)
    s_in = lastindex(input)
    
    for n=eachindex(fir)
        So_conv[n:s_in+n-1] .+= fir[n] .* input
    end

    return nothing
end

function u_gen_ir_rc(dt,bw,t_len)
    tt = [0:dt:t_len-dt;]

    ω = (2*π*bw)
    ir = ω*exp.(-tt*ω)
    ir = ir/sum(ir*dt)

    return ir
end

function u_hist(samples, edges)
    h = fit(Histogram, samples, edges)
    return h.weights
end

function u_hist(samples, minval, maxval, nbin)
    weights = zeros(nbin)
    bin_size = (maxval-minval)/nbin

    for s in samples
        idx = Int(floor((s-minval)/bin_size))+1
        idx = idx < 1 ? 1 : idx > nbin ? nbin : idx
        weights[idx] += 1
    end
    return weights
end

function u_fr_to_imp(f, H, Tsym, osr; npre=50, npost=200, savename = "", t_name="dt", ir_name="ir")
    fbaud = 1/Tsym
    fbaud_max = 2*maximum(f)

    @assert fbaud < fbaud_max "Max frequency too low for desired symbol rate"

    fbaud_max = fbaud*floor(fbaud_max/fbaud)

    Hm = abs.(H)
    Hp = unwrap(angle.(H))

    if f[1] == 0
        Hm_ds = [reverse(Hm[2:end-1]); Hm]
        Hp_ds = [reverse(-Hp[2:end-1]; Hp)]
        f_ds = [-reverse(f[2:end-1]); f]
    else
        Hm_ds = [reverse(Hm[1:end-1]); abs(Hm[1]); Hm]
        Hp_ds = [reverse(-Hp[1:end-1]); 0; Hp]
        f_ds = [-reverse(f[1:end-1]); 0; f]
    end
    
    num_fft_pts = 2^16
    df = fbaud_max/2/num_fft_pts

    f_ds_itp = -fbaud_max/2+df:df:fbaud_max/2
    itp_Hm = linear_interpolation(f_ds, Hm_ds)
    Hm_ds_itp = ifftshift(itp_Hm.(f_ds_itp))
    itp_Hp = linear_interpolation(f_ds, Hp_ds)
    Hp_ds_itp = ifftshift(itp_Hp.(f_ds_itp))

    ir = fbaud_max*real.(ifft(Hm_ds_itp.*exp.(im.*Hp_ds_itp)))

    dt = 1/fbaud_max
    tt = 0:dt:(length(ir)-1)*dt
    itp_ir = linear_interpolation(tt,ir);
    
    dt_s = Tsym/osr
    tt_s =0:Tsym/osr:tt[end]
    ir_itp = itp_ir.(tt_s)
    
    max_idx = argmax(ir_itp)
    start_idx = (max_idx-npre*osr > 0) ? max_idx-npre*osr : 1
    end_idx = (max_idx+npost*osr-1 < length(ir_itp)) ? max_idx+npost*osr-1 : length(ir_itp)

    ir_itp = ir_itp[start_idx:end_idx]

    if savename != ""
        matwrite(savename, Dict(
            dt_name => dt_s,
            ir_name => ir_itp
        ))
    end

    return ir_itp
end

function u_fr_to_imp(filename, Tsym, osr; npre=50, npost=200, savename = "", freq_name="fr", tf_name="H")
    data = matopen(filename)
    fr = vec(read(data, freq_name))
    H = vec(read(data, tf_name))
    close(data)
    return u_fr_to_imp(fr, H, Tsym, osr, npre=npre, npost=npost, savename = savename)
end

end