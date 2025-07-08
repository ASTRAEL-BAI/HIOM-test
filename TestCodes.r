g = make_lattice(c(5,5,1),nei =1)

plot(g, layout = layout_on_grid(g))
windows()
distances(g)