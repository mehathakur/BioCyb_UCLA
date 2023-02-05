function plot_sim(sol; title::AbstractString = "Thyrosim simulation", automargins::Bool=true)

    # parameters to adjust figure limits
    p = sol.prob.p 
    t4lim, t3lim, tshlim = 140, 4, 10
    T4 = 777.0 * sol[1, :] / p[47]
    T3 = 651.0 * sol[4, :] / p[47]
    TSH = 5.6 * sol[7, :] / p[48]
    if automargins
        t4lim = max(1.2maximum(T4), 130.0)
        t3lim = max(1.2maximum(T3), 2.5)
        tshlim = max(1.2maximum(TSH), 5.5)
    end

    p1 = plot(sol.t / 24.0, T4, ylim=(0, t4lim), label="",
       ylabel="T4 (mcg/L)", title=title)
    p1 = hline!([45, 120], label= "")
    
    p2 = plot(sol.t / 24.0, T3, ylim=(0, t3lim), label="", 
       ylabel="T3 (mcg/L)")
    p2 = hline!([0.6, 1.8], label= "")
    
    p3 = plot(sol.t / 24.0, TSH, ylim=(0, tshlim), label="",
       ylabel="TSH (mU/L)", xlabel="time [days]")
    p3 = hline!([0.45, 4.5], label= "")
    
    plot(p1, p2, p3, layout=(3, 1))
end