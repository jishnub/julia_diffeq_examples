using DifferentialEquations
using DiffEqPhysics
using PyPlot

function H(p,q,params)
	m = params[1]
	l = params[2]
	g = params[3]

	p^2/(2m*l^2) + m*g*l*(1-cos(q))
end

function main()
	m=1.;l=1.;g=9.8;

	params = [m,l,g]

	tspan = (0.0,10.0)

	p0 = 0.
	q0 = π/6

	prob = HamiltonianProblem(H,p0,q0,tspan,params)

	sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

	plot(sol.t,sol[1,1:end],label="velocity")
	plot(sol.t,sol[2,1:end],label="displacement")

	axhline(0,ls="dotted",color="grey",zorder=0)
	axhline(q0,ls="dotted",color="grey",zorder=0)
	axhline(-q0,ls="dotted",color="grey",zorder=0)

	legend(loc="best")
	margins(y=0.2)

end

main()