################################################################################

# using  ApproXD, Dierckx, Grid, NumericalMath, PyPlot, PyCall
using  ApproXD, Grid

################################################################################
# Interpolation schemes to consider:
methods_1D = ["Grid (linear)", "ApproXD (linear)"]
methods_2D = ["Grid (linear)", "ApproXD (linear)"]

# Function to interpolate
f(x::Float64) = -(x^(-1.0))

# Grids for interpolation (regular/irregular)
zmin = -0.1
zmax = 10.0
ymin = 0.2
ymax = 30.0
xpoints = 100
ypoints = 20
zpoints = 15

xmin = ymin + zmin
xmax = ymax + zmax

xgrid = linspace(xmin, xmax, xpoints)
ygrid = linspace(ymin, ymax, ypoints)
zgrid = linspace(zmin, zmax, zpoints)

################################################################################
# Evaluate functions on grid points

V_1D = Array(Float64, xpoints)
V_2D = Array(Float64, (ypoints, zpoints))

# Values on grid, 1 dimension
for i = 1:xpoints
  V_1D[i] = f(xgrid[i])
end

# Values on grid, 2 dimensions
for i = 1:ypoints
  for j = 1:zpoints
    V_2D[i, j] = f(ygrid[i] + zgrid[j])
  end
end

################################################################################
# Interpolate in 1 dimension

# Grid (regular grid only)
xrange = range(xgrid[1], xgrid[2]-xgrid[1], xpoints)

V_1D_int_Lin =
  CoordInterpGrid(xrange, V_1D, BCnearest, InterpLinear)

# ApproXD (regular/irregular grid)
xs = Array{Float64, 1}[]
push!(xs, xgrid)
V_1D_int_ApproXD = Lininterp(V_1D, xs)

################################################################################
# Interpolate in 2 dimensions

# Grid
yrange = range(ygrid[1], ygrid[2]-ygrid[1], ypoints)
zrange = range(zgrid[1], zgrid[2]-zgrid[1], zpoints)

V_2D_int_Lin =
  CoordInterpGrid((yrange, zrange), V_2D, BCnearest, InterpLinear)

# ApproXD (regular/irregular grid)
yz = Array{Float64, 1}[]
push!(yz, ygrid)
push!(yz, zgrid)
V_2D_int_ApproXD = Lininterp(V_2D, yz)

################################################################################
# Evaluate Regular Grid Interpolants Off-Grid

xop = 200*xpoints  # Number of off-grid points
yop = 200*ypoints
zop = 200*zpoints

xoffgrid = linspace(xgrid[1], xgrid[end], xop)
yoffgrid = linspace(ygrid[1], ygrid[end], yop)
zoffgrid = linspace(zgrid[1], zgrid[end], zop)

V_1D_offgrid = Array(Float64, (length(methods_1D), xop))
V_1D_actual = similar(V_1D_offgrid)

V_2D_offgrid = Array(Float64, (length(methods_2D), yop, zop))
V_2D_actual = similar(V_2D_offgrid)


g1d = @elapsed for i = 1:xop
  # Grid
  V_1D_offgrid[1, i] = V_1D_int_Lin[xoffgrid[i]]
end

a1d = @elapsed for i = 1:xop
  # ApproXD (regular/irregular grid)
  V_1D_offgrid[2, i] = getValue(V_1D_int_ApproXD, [xoffgrid[i]])[1]
end


g2d = @elapsed for i = 1:yop
  for j = 1:zop
    yt = yoffgrid[i]
    zt = zoffgrid[j]
    # Grid
    V_2D_offgrid[1, i, j] = V_2D_int_Lin[yt, zt]
  end
end

a2d = @elapsed for i = 1:yop
  for j = 1:zop
    yt = yoffgrid[i]
    zt = zoffgrid[j]
    # ApproXD (regular/irregular grid)
    V_2D_offgrid[2, i, j] = getValue(V_2D_int_ApproXD, [yt, zt])[1]
  end
end


println("timing for 1D Grid: $g1d")
println("timing for 1D ApproXD: $a1d")
println("timing for 1D ApproXD: $g2d")
println("timing for 1D ApproXD: $a2d")

# ######################################################################

# # Calculate Errors
# Error_1D = (V_1D_actual - V_1D_offgrid)./V_1D_actual
# Error_2D = (V_2D_actual - V_2D_offgrid)./V_2D_actual

# ######################################################################
# # Print results
# @printf "\n∑|f(x)-f(x_int)/f(x)|/n in 1 dimension\n"
# for i = 1:length(methods_1D)
#   @printf "\t%s: %.4f\n" methods_1D[i] sum(abs(Error_1D[i, :]))/xop
# end

# @printf "\n∑|f(x)-f(x_int)/f(x)|/n in 2 dimensions\n"
# for i = 1:length(methods_2D)
#   @printf "\t%s: %.4f\n" methods_2D[i] sum(abs(Error_2D[i, :, :]))/(yop*zop)
# end

# # Plot 1D results
# fig, ax = PyPlot.subplots(2,1, figsize = (12,10))
# for i = 1:7
#   ax[1,1][:plot](V_1D_offgrid[i, 1:end/20]', label = methods_1D[i])
#   ax[2,1][:plot](Error_1D[i, 1:end/20]')
# end
# ax[1,1][:plot](V_1D_actual[1, 1:end/20]', label = "True Value",
#                linestyle = ":", color = "black")
# ax[1,1][:set_title]("True Value vs. Interpolated Values")
# ax[1,1][:legend](loc = "best")
# ax[2,1][:set_title]("Relative Error")
# ax[2,1][:legend](loc = "best")
# fig[:suptitle]("Interpolation in One Dimension", fontsize = 16)
# plt.show()

