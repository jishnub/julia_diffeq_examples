using DifferentialEquations
using PyPlot

function rhs(u::Array{Float64,3},p::Vector{Float64},t::Float64)

	disp = u[:,:,1]
	vel = u[:,:,2]

	dx = p[1]
	dy = p[2]
	c = p[3]

	dtvel = c^2*∇2(disp,dx,dy)  ::Array{Float64,2}

	du = cat(3,vel,dtvel)

end

function ddx(f::Array{Float64,2},dx::Float64,order::Int64=1)
	Nx = size(f,2)
	Lx = Nx*dx
	Nk = div(Nx,2)+1
	kx = 2π/Lx*collect(0:Nk-1)
	df = irfft(((im*kx).^order)'.*rfft(f,2),Nx,2)
end

d2dx2(f::Array{Float64,2},dx::Float64) = ddx(f,dx,2)

function ddy(f::Array{Float64,2},dy::Float64,order::Int64=1)
	Ny = size(f,1)
	Ly = Ny*dy
	Nk = div(Ny,2)+1
	ky = 2π/Ly*collect(0:Nk-1)
	df = irfft(((im*ky).^order).*rfft(f,1),Ny,1)
end

d2dy2(f::Array{Float64,2},dx::Float64) = ddy(f,dx,2)

∇2(f::Array{Float64,2},dx::Float64,dy::Float64) = d2dx2(f,dx) + d2dy2(f,dy)

function main()
	
	σx = 1.; x0=4.; kx = 2π; ky = 0.; c = 1.; ω=c*kx; Lx=10.; Ly=10.
	x = linspace(0,Lx,65)[1:end-1] # assume periodicity, leave out last gridpoint
	y = linspace(0,Ly,65)[1:end-1]
	dx = x.step.hi # step size
	dy = y.step.hi

	disp0 = ones(y).*(exp.(-(x-x0).^2/2σx^2).*sin.(kx*(x-x0)))';
	vel0 = ones(y).*(exp.(-(x-x0).^2/2σx^2).*(-ω).*cos.(kx*(x-x0)))';

	u0 = cat(3,disp0,vel0)

	tspan = (0.0,5.0)

	prob = ODEProblem(rhs,u0,tspan,[dx,dy,c])

	sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

	# plot the wave displacement
	if !isdir("plots")
		mkdir("plots")
	end

	for timeind in 1:10:length(sol.t)
		pcolormesh(x,y,sol[timeind][:,:,1])
		savefig("plots/$timeind.png")
		cla()
	end
end

main()