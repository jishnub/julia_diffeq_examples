using DifferentialEquations
using PyPlot

function rhs(u::Array{Float64,2},p::Vector{Float64},t::Float64)

	disp = u[:,1] ::Vector{Float64}
	vel = u[:,2] ::Vector{Float64}

	dx = p[1] ::Float64
	c = p[2] ::Float64

	dtvel = c^2*d2dx2(disp,dx)  ::Vector{Float64}

	du = hcat(vel,dtvel)

end

function ddx(f::Vector{Float64},dx::Float64,order::Int64=1)
	Nx = length(f)
	Lx = Nx*dx
	Nk = div(Nx,2)+1
	kx = 2π/Lx*collect(0:Nk-1)
	df = irfft(((im*kx).^order).*rfft(f),Nx)
end

d2dx2(f::Vector{Float64},dx::Float64) = ddx(f,dx,2)

function main()
	
	σx = 1; x0=4; kx = 2π/2; c = 1; ω=c*kx; Lx=10
	x = linspace(0,Lx,129)[1:end-1] # assume periodicity, leave out last gridpoint
	dx = x.step.hi # step size

	disp0 = exp.(-(x-x0).^2/2σx^2).*sin.(kx*(x-x0));
	vel0 = exp.(-(x-x0).^2/2σx^2).*(-ω).*cos.(kx*(x-x0));

	u0 = hcat(disp0,vel0)

	tspan = (0.0,5.0)

	prob = ODEProblem(rhs,u0,tspan,[dx,c])

	sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

	# plot the wave displacement
	if !isdir("plots")
		mkdir("plots")
	end

	for timeind in 1:10:length(sol.t)
		plot(x,sol[timeind][:,1])
		savefig("plots/$timeind.png")
		cla()
	end
end

main()