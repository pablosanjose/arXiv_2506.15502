# using GLMakie
using CairoMakie
using MakieTeX
using ProgressMeter
using ArnoldiMethod, Arpack, LinearMaps
using JLD2
using Quantica
using ColorSchemes
using Revise
includet("../code.jl")

## Figure 2
    ## Data
        function data_figure2(nev = 8, npoints = 30; params...)
            (; hN, ϵ1) = build(; params...);
            _, psi = lowest_eigenpair(hN(; params...), nev)
            ns = [range(0, 0.49999, npoints+1); range(0.5000001, 1.49999, 2npoints + 1)]
            solver = ES.ShiftInvert(ES.Arpack(nev = 4+nev), 0)
            b = bands(hN, ns; mapping = n -> ftuple(; params..., n), solver)
            Δ = NamedTuple(params).Δ
            params´ = (; params..., Δ0 = Δ, Δ1 = 0.9*Δ)
            b´ = bands(hN, range(-1,1,2npoints+1); mapping = δΔ -> ftuple(; params´..., Δ0 = Δ - Δ*δΔ/2, Δ1 = Δ + Δ*δΔ/2), solver)
            nt = NamedTuple(params);
            return (; psi, ns, b, b´, ϵ1, params = nt)
        end

        function lowest_eigenpair(h, nev = 8)
            es, psis = spectrum(h, solver = ES.ShiftInvert(ES.ArnoldiMethod(nev = 4+nev), 0))[1:nev, around = 0]
            i = findfirst(<(0), real(es))
            psi = reinterpret(SVector{4,ComplexF64}, psis[:,i])
            return es[i], psi
        end

        params = (; R = 1000, µ = 2, ϕ = 0, n = 1, Δ = 0.2, α = 0, Vz = 0.0001)
        data1 = data_figure2(64; params...)

        params = (; R = 1000, µ = 4, ϕ = 0, n = 1, Δ = 0.2, α = 0, Vz = 0.0001)
        data2 = data_figure2(72; params...)

        params = (; R = 1000, µ = 4, ϕ = 0, n = 0.5000001, Δ = 0.2, α = 0, Vz = 0.0001)
        data3 = data_figure2(72; params...)

        jldsave("figures/fig2_data.jld2"; data1, data2, data3)

    ## Plot
        @load "figures/fig2_data.jld2" data1 data2

        function figure2(data, red´ = 37:40)
            (; psi, b, b´, ϵ1, params) = data
            Δ = params.Δ
            common = (; xlabelsize = 16, ylabelsize = 16)
            labelcommon = (; fontsize = 16,
                padding = (0, 20, -10, 0),
                halign = :right)

            ## (a) ##
            img = PDFDocument(read("figures/Fig2a.pdf"));
            fig = Figure(size = (600,550))
            ax = Axis(fig[1,1]; alignmode = Outside(-60,-60,-60,-20))
            teximg!(ax, img; scale = 1.5)
            hidedecorations!(ax)
            hidespines!(ax)
            Label(fig.layout[1, 1, TopLeft()], "(a)"; labelcommon...)

            ## (b) ##
            ϕs = range(-1, 1, length(psi))
            ax = Axis(fig[2,1];
                xlabel = L"\varphi", xticks = (-1:0.5:1, ["-π","-π/2","0","π/2","π"]),
                ylabel = L"\text{wavefunction (arb. units)}", yticks = (0:0,["0"]),
                limits = ((-1,1), nothing), ygridvisible = false, common...)
            lines!(ax, ϕs, norm.(psi); color = :blue, label = L"|\Psi|^2")
            lines!(ax, ϕs, real.((x->x[4]).(psi)); color = :orange, linewidth = 0.6, label = L"\text{Re}\Psi_e")
            lines!(ax, ϕs, real.((x->x[1]).(psi)); color = :purple, linewidth = 0.6, label = L"\text{Re}\Psi_h")
            scatter!(0.9, 0.068; color = :blue)
            axislegend(ax; position = :lt)
            Label(fig.layout[2, 1, TopLeft()], "(b)"; labelcommon...)

            ## (c) ##
            ax = Axis(fig[1,2];
                title = L"\Delta_0=\Delta_1", titlealign = :right,
                xlabel = L"\Phi/\Phi_0", ylabel = L"\epsilon/(\Delta_0+\Delta_1)",
                yticks = (range(-2Δ,2Δ,3),["-1","0","1"]),
                limits = ((0, 1.5), (-2.5Δ, 2.5Δ)), xgridvisible = false, common...)
            cs = ColorScheme([colorant"red", colorant"grey"])
            qplot!(b; color = (ψ, ϵ, ϕs) -> ifelse(abs(ϵ)<0.12,0,1), colormap = cs, hide = :nodes)
            scatter!(1, ϵ1; color = :blue)
            text!(ax, 1.3, -0.11; text = L"\times 2", color = :red, fontsize = 12)
            text!(ax, 1.3, 0.05; text = L"\times 2", color = :red, fontsize = 12)
            vlines!(ax, [0.5]; color = :gray, linewidth = 0.5)
            Label(fig.layout[1, 2, TopLeft()], "(c)"; labelcommon...)

            ## (d) ##
            ax = Axis(fig[2,2];
                title = L"\Phi/\Phi_0 = 1", titlealign = :right,
                xlabel = L"(\Delta_1-\Delta_0)/\Delta_0", ylabel = L"\epsilon/\Delta_0",
                yticks = (range(-2Δ,2Δ,3),["-0.5","0","0.5"]),
                limits = ((-1, 1), (-2.5Δ, 2.5Δ)))
            qplot!(b´; color = :grey, hide = :nodes)
            qplot!(b´[red´]; color = :red, hide = :nodes)
            scatter!(0, ϵ1; color = :blue)
            text!(ax, 0.7, -0.16; text = L"\times 2", color = :red, fontsize = 12)
            text!(ax, 0.7, 0.1; text = L"\times 2", color = :red, fontsize = 12)
            Label(fig.layout[2, 2, TopLeft()], "(d)"; labelcommon...)

            colsize!(fig.layout, 2, Auto(0.6))
            rowgap!(fig.layout, 5)
            colgap!(fig.layout, 25)
            return fig
        end

        f = figure2(data2, 37:40)
        save("figure2.pdf", f)

## Figure 3
    ## Data
        function data_figure3(αs, ns = range(0.6, 1.4, 9); params...)
            (; hN, p) = build(; params...);
            nt = NamedTuple(params)
            ps = @showprogress [[operator(hN, p; params..., α, n) for α in αs] for n in ns]
            return (; ps, αs, ns, R = nt.R)
        end

        params = (; R = 500, µ = 4, ϕ = 0, n = 0.9, Δ = 0.2, Vz = 0.001)
        data1 = data_figure3(range(0, 5, 151); params...)

        params = (; R = 750, µ = 4, ϕ = 0, n = 0.9, Δ = 0.2, Vz = 0.001)
        data2 = data_figure3(range(0, 5, 151); params...)

        params = (; R = 1000, µ = 4, ϕ = 0, n = 0.9, Δ = 0.2, Vz = 0.001)
        data3 = data_figure3(range(0, 5, 151); params...)

        params = (; R = 1500, µ = 4, ϕ = 0, n = 0.9, Δ = 0.2, Vz = 0.001)
        data4 = data_figure3(range(0, 5, 151); params...)

        # Data SU(2)
        function λdata(idat, data, n, α)
            δα = data.αs[2]-data.αs[1]
            δn = data.ns[2]-data.ns[1]
            i = findfirst(x -> n - δn/2 <= x <= n+δn/2, data.ns)
            j = findfirst(x -> α - δα/2 <= x <= α+δα/2, data.αs)
            pmat = data.ps[i][j]
            R = data.R
            # (1/2)*2πRλ⁻¹ᵢ = pᵢ * R (Eq. 6), pᵢ = 1/2*Tr(p*σᵢ)
            λp⁻¹2πR = 2R*SU2toSO3(pmat)[2]
            λz⁻¹2πR = 2R*SU2toSO3(pmat)[4]
            return idat, (data.ns[i], data.αs[j]), (λp⁻¹2πR, λz⁻¹2πR)
        end

        λs = λdata(4, data4, 1.2, 2.25), λdata(4, data4, 0.9, 1.5), λdata(3, data3, 1.0, 3.0)

        using OrdinaryDiffEq

        # pp, pz = λp⁻¹2πR, λz⁻¹2πR,
        function dataode_figure3((pp, pz), ϕmax = 2π)
            rot(u, p, ϕ) = -0.5 * im * SA[pz pp*cis(-ϕ); pp*cis(ϕ) -pz] * u
            prob = ODEProblem(rot, SA[1.0 + 0im, 0], (0.0, ϕmax))
            sol = solve(prob, Tsit5())
            return sol
        end

        dataode1 = dataode_figure3(λs[1][3], 40π)
        dataode2 = dataode_figure3(λs[2][3], 40π)
        dataode3 = dataode_figure3(λs[3][3], 40π)

        jldsave("figures/fig3_data.jld2"; data1, data2, data3, data4, λs)

    ## Plots
        @load "figures/fig3_data.jld2" data1 data2 data3 data4 λs # dataode1 dataode2 dataode3

        function figure3(datas, dataodes, λs)
            fig = Figure(size = (600, 700))

            # (p,z) := (1/2)*2πRλ⁻¹ᵢ = pᵢ * R (Eq. 6), pᵢ = 1/2*Tr(p*σᵢ)
            pRz(pmat) = SU2toSO3(pmat)[4]
            pRp(pmat) = SU2toSO3(pmat)[2]

            common = (; xlabelsize = 16, ylabelsize = 16)
            labelcommon = (; fontsize = 16, halign = :right, padding = (0, 30, -10, 0))
            orbitcolors = (:gray, :orange, Makie.RGBA(0.7,0,0.7))
            # orbitcolors = (:gray, :orange, Makie.RGBA(0.3, 0.7, 0))

            ## (a-d) ##
            for (n, (r, c), data, tag) in zip(eachindex(datas), ((1,1), (1,2), (2,1), (2,2)), datas, ("(a)", "(b)", "(c)", "(d)"))
                ax = Axis(fig[r, c];
                    xlabel = L"\alpha \text{[meV⋅nm]}", ylabel = L"2\pi R \lambda^{-1}_i",
                    title = "R = $(data.R) nm",
                    # limits = (nothing, (-3,2.5)),
                    common...)
                for row in data.ps
                    lines!(ax, data.αs, pRp.(row) .* (2*data.R); color = :blue, linewidth = 1, label = L"2\pi R\lambda^{-1}_\parallel")
                end
                for row in data.ps
                    lines!(ax, data.αs, pRz.(row) .* (2*data.R); color = :red, linewidth = 1, label = L"2\pi R\lambda^{-1}_z")
                end
                hlines!(ax, 1; color = :red, linestyle = :dash, linewidth = 1)
                for ((i, (_, α), (pp, pz)), color) in zip(λs, orbitcolors)
                    if n == i
                        scatter!(α, pz; color = :black, markersize = 10)
                        scatter!(α, pz; color, markersize = 8)
                        scatter!(α, pp; color = :black, markersize = 10)
                        scatter!(α, pp; color, markersize = 8)
                    end
                end
                padding = (n == 2 || n == 4) ? (0, 10, -10, 0) : (0, 30, -10, 0)
                Label(fig.layout[r, c, TopLeft()], tag; labelcommon..., padding)
                # n == 2 && text!(0, 0.5; text = "optimal", color = :red)
                n == 3 && axislegend(ax; position = :lb, unique = true)
                n == 4 && ylims!(ax, -5.4, 5.4)
                n == 4 && text!(0.6, 0.0; text = "optimal", color = :red)
                n == 2 && (arrows!(ax, [1], [-5], [0], [-4]); text!(ax, 0.1, -9; text = L"\Phi/\Phi_0"))
                c == 2 && hideydecorations!(ax; grid = false, ticks = false, ticklabels = false)
                r == 1 && hidexdecorations!(ax; grid = false, ticks = false)
            end

            gl = fig[3,1:2] = GridLayout()
            ax = Axis(gl[1,1];
                xlabel = L"\varphi_B", ylabel = L"\cos\,\theta_B",
                limits = ((-π,π), (-1.04, 1.04)),
                xticks = (range(-π,π,5), ["-π","-π/2","0","π/2","π"]),
                common...)
            for (dataode, color) in zip(dataodes, orbitcolors)
                φs = range(extrema(dataode.t)..., 3000)
                zs, ϕs = Float64[],Float64[]
                ϕprev = 0.0
                zprev = 1.0
                for φ in φs
                    ϕ, z = SU2toPolar(dataode(φ))
                    if abs(ϕ-ϕprev) > 0.9 * 2π
                        push!(zs, (zprev+z)/2, NaN, (zprev+z)/2)
                        push!(ϕs, -sign(ϕ) * π, NaN, sign(ϕ) * π)
                    end
                    ϕprev = ϕ
                    zprev = z
                    push!(zs, z)
                    push!(ϕs, ϕ)
                end
                lines!(ax, ϕs, zs; color)
            end
            Label(gl[1, 1, TopLeft()], "(e)"; labelcommon...)

            ax = Axis3(gl[1,2]; aspect = :data, alignmode = Outside(-60,-60,-60,-60))
            for (dataode, color) in zip(dataodes, orbitcolors)
                φs = range(extrema(dataode.t)..., 1500)
                xs, ys, zs = ([SU2toSO3(dataode(φ))[i] for φ in φs] for i in 1:3)
                lines!(xs, ys, zs; color)
            end
            hidedecorations!(ax)
            hidespines!(ax)
            text!(-1,0,-1.3; text = "Bloch sphere")
            colsize!(gl, 2, Auto(0.3))

            return fig
        end

        f = figure3((data1, data2, data3, data4), (dataode1, dataode2, dataode3), λs)
        save("figure3.pdf", f)

## Figure 4
    ## Data
        function data_figure4(xs = range(-1.1,1.1,100); nev = 8, params...)
            (; Δ0, Δ1, α, R, n) = NamedTuple(params)
            (; hN, t) = build(; params...)
            ## bands ##
            solver = ES.ShiftInvert(ES.ArnoldiMethod(; nev), 0)
            b = bands(hN, xs;
                mapping = x -> ftuple(; params..., µN = 2t*(1+x)),
                split = false, solver)
            ## eigenvalues
            solver = ES.ShiftInvert(ES.ArnoldiMethod(; nev = 4), 0)
            ϵs = @showprogress [first(pick_subspace(
                spectrum(hN; solver, params..., µN = 2t*(1-x)))) for x in xs]

            return (; xs, b, ϵs, α, R, n, t, Δ0, Δ1)
        end

        params = (; R = 1000, N = 1000, ϕ = 0, n = 0.50001, Δ1 = 0.2, Δ0 = 0.2, α = 3, Vz=0.0)
        data1 = data_figure4(range(-1.01,1.01,400); nev = 48, params...)

        params = (; R = 1000, N = 1000, ϕ = 0, n = 0.505, Δ1 = 0.2, Δ0 = 0.2, α = 3, Vz=0.0)
        data2 = data_figure4(range(-1.01,1.01,400); nev = 48, params...)

        params = (; R = 1000, N = 1000, ϕ = 0, n = 0.55, Δ1 = 0.2, Δ0 = 0.2, α = 3, Vz=0.0)
        data3 = data_figure4(range(-1.01,1.01,400); nev = 48, params...)

        params = (; R = 1000, N = 1000, ϕ = 0, n = 0.85, Δ1 = 0.2, Δ0 = 0.25, α = 3, Vz=0.0)
        data4 = data_figure4(range(-1.01,1.01,400); nev = 48, params...)

        jldsave("figures/fig4_data.jld2"; data1, data2, data3, data4)

    ## Plots
        @load "figures/fig4_data.jld2" data1 data2 data3 data4

        function figure4(data1, data2, data3, data4)
            fig = Figure(size = (600, 500))

            common = (; xlabelsize = 16, ylabelsize = 16, xgridvisible = false)
            labelcommon = (; fontsize = 16, padding = (0, 40, 0, 0), halign = :right)

            ## (a) ##
            ax = Axis(fig[1, 1];
                title = L"\Delta_0 = \Delta_1,\,\, \Phi/\Phi_0 \approx 0.5",
                xlabel = L"\mu/W", ylabel = L"\epsilon \text{[meV]}",
                limits = (nothing, (-0.5,0.5)), common...)
            cs = ColorScheme([colorant"red", colorant"gray"])
            qplot!(data1.b; color = (ψ, ϵ, ϕs) -> ifelse(abs(ϵ)<0.05,0,1), colormap = cs, hide = :bands, nodesizefactor = 2)
            text!(ax, 0.7, 0.01; text = L"\times 4", color = :red, fontsize = 12)
            Label(fig.layout[1, 1, TopLeft()], "(a)"; labelcommon...)

            ## (b) ##
            ax = Axis(fig[1, 2];
                title = L"\Delta_0 > \Delta_1,\,\, \Phi/\Phi_0 > 0.5",
                xlabel = L"\mu/W", ylabel = L"\epsilon \text{[meV]}",
                limits = (nothing, (-0.5,0.5)), common...)
            cs = ColorScheme([colorant"red", colorant"gray"])
            qplot!(data4.b; color = (ψ, ϵ, ϕs) -> ifelse(abs(ϵ)<0.05,0,1), colormap = cs, hide = :bands, nodesizefactor = 2)
            text!(ax, 0.7, 0.04; text = L"\times 2", color = :red, fontsize = 12)
            text!(ax, 0.7, -0.09; text = L"\times 2", color = :red, fontsize = 12)
            hideydecorations!(ax, ticks = false)
            Label(fig.layout[1, 2, TopLeft()], "(b)"; labelcommon..., padding = (0, 10, 0, 0))

            ## (c) ##
            ax = Axis(fig[2, 1:2];
                title = L"\text{residual splitting at $\Delta_0 = \Delta_1$}",
                xlabel = L"\mu/W", ylabel = L"\delta\epsilon^\alpha \text{[$\mu$eV]}",
                # limits = (nothing, (-0.05,-0.022)),
                ygridvisible = false,
                common...)
            colors = (:blue, :green, :orange)
            for (data, color) in zip((data1, data2, data3), colors)
                (; xs, ϵs, n) = data
                δϵ = (((ϵ1,ϵ2),) -> 1e3*abs(ϵ1-ϵ2)).(ϵs)
                lines!(ax, xs, δϵ; color, label = "$(round(n, digits = 3))")
            end
             (; xs, ϵs, α, R, n) = data3
            δϵa = [1e3*α/R * n * abs(x) for x in xs]
            lines!(ax, xs, δϵa, color = :black, linestyle = :dash)
            axislegend(ax, L"\Phi/\Phi_0"; position = :ct, orientation = :horizontal)
            Label(fig.layout[2, 1, TopLeft()], "(c)"; labelcommon...)

            rowsize!(fig.layout, 2, Auto(0.5))
            return fig
        end

        f = figure4(data1, data2, data3, data4)
        save("figure4.pdf", f)
