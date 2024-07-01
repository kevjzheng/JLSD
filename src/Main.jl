using Revise, BenchmarkTools, ProfileCanvas
# import Plots as plt; pythonplot()

include("./tb/TB.jl")


trx = init_trx()
# VSCodeServer.@profview_allocs run_sim(trx, run_sim_blk)
@time run_blk_iter(trx, 0, trx.param.nblk, sim_blk)

trx.wvfm.en_plot=true
w_plot_test(trx.wvfm)
update_eye(trx.wvfm.eye1, x_nui=2, x_npts_ui=256, y_npts=256, x_ofst_ui = 0)#trx.cdr.pi_code/trx.clkgen.pi_codes_per_ui)
reset_limits!(trx.wvfm.axes[3,2])


ber = trx.bist.ber_err_cnt/trx.bist.ber_bit_cnt
println("BER = $ber")
    
# a = trx.drv.buffer_debug[10000:end];
# b = u_unwrap_0x(a);

# f1 = w_newfig()
# # lines!(Axis(f1[1,1]),a)
# lines!(Axis(f1[1,1]),b)
# jitter_bnd = (-0.1,0.1).+(-6, 6).*(trx.drv.rj_s/trx.param.tui) .+ (-1.2, 1.2).*trx.drv.sj_amp_ui
# density!(Axis(f1[1,2]), b .- mean(b), boundary=jitter_bnd, npoints=200)


 
