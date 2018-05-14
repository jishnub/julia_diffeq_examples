using DifferentialEquations
using DiffEqPhysics
using Plots
# import PyPlot
# plt = PyPlot

function H(p,q,params)
	m1 = params[1]
	m2 = params[2]
	l1 = params[3]
	l2 = params[4]
	g = params[5]

	(p[1]^2 * l2^2*m2 + l1^2*(m1+m2)*p[2]^2 - 2m2*l1*l2*p[1]*p[2]*cos(q[1]-q[2]))/
	(2l1^2*l2^2*m2*(m1+m2*sin(q[1]-q[2])^2)) - 
	m2*g*l2*cos(q[2]) - (m1+m2)*g*l1*cos(q[1])
end

function main()

	m1=1.;m2=2.;l1=2.;l2=1.5;g=9.8;
	params=[m1,m2,l1,l2,g]

	tspan = (0.,10.)

	p1_0 = 0.
	q1_0 = π/6

	p2_0 = 0.
	q2_0 = π/3

	p0 = [p1_0,p2_0]
	q0 = [q1_0,q2_0]

	prob = HamiltonianProblem(H,p0,q0,tspan,params)

	sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
	
	anim = @animate for timeind = 1:2:length(sol.t)

		height = 2*(l1+l2)
		θ1 = sol[3,timeind]
		θ2 = sol[4,timeind]
		x1 = l1*sin(θ1); x2 = l1*sin(θ1) + l2*sin(θ2);
		y1 = height - l1*cos(θ1); y2 = height - (l1*cos(θ1) + l2*cos(θ2));

		p = plot(xlims=(-1.2*(l1+l2),1.2*(l1+l2)),ylims=(0,1.2*height),legend=false)
		plot!([0],[height],seriestype=:scatter,marker=3)
		plot!([x1],[y1],seriestype=:scatter)
		plot!([x2],[y2],seriestype=:scatter)
		plot!([0,x1],[height,y1])
		plot!([x1,x2],[y1,y2])

	end

	gif(anim,"plots/double_pendulum.gif",fps=12)


end

main()