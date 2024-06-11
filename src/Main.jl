using Revise, BenchmarkTools, ProfileCanvas
# import Plots as plt; pythonplot()

include("./tb/TB.jl")


trx = init_trx()
@time run_sim(trx, run_sim_blk)

trx.wvfm.en_plot=true
w_plot_test(trx.wvfm)
update_eye(trx.wvfm.eye1, x_nui=1, x_npts_ui=256, y_npts=256, x_ofst_ui = trx.cdr.pi_code/trx.clkgen.pi_codes_per_ui)
reset_limits!(trx.wvfm.axes[3,2])


ber = trx.bist.ber_err_cnt/trx.bist.ber_bit_cnt
println("BER = $ber")
    
# a = trx.drv.buffer_debug[1000:end];
# b = u_unwrap_0x(a);
# c = u_hist(b.-mean(b), -0.5:0.001:0.5)

# f1 = w_newfig()
# lines!(Axis(f1[1,1]),a)
# lines!(Axis(f1[1,2]),b)
# density!(Axis(f1[1,3]),b .- mean(b))


