ip_geo = Lagrange{1,RefCube,1}()
ip_f = Lagrange{1,RefCube,2}()
ip_jump = FerriteCon2022.JumpInterpolation(ip_f)
ip_mid = FerriteCon2022.MidPlaneInterpolation(ip_geo)

qr = QuadratureRule{1,RefCube}(2)

cv = FerriteCon2022.CohesiveVectorValues(qr, ip_jump, ip_mid)

xe = Vec.([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)])
reinit!(cv, xe)

xs = Vec.([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (2.0, 0.0), (2.0, 1.0)])


nodes = Node.(xs)
cells = [FerriteCon2022.CohesiveQuadrilateral((1,2,3,4)), FerriteCon2022.CohesiveQuadrilateral((3,5,4,6))]
grid = Grid(cells, nodes)

dh = MixedDofHandler(grid)
push!(dh, :u, 2, ip_jump)
close!(dh)

ip_base = Lagrange{2, RefCube, 1}()
ip_jump = FerriteCon2022.JumpInterpolation(ip_base)
qr = QuadratureRule{2,RefCube}(2)
cv = FerriteCon2022.CohesiveVectorValues(qr, ip_jump)

xe = Vec.([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)])
reinit!(cv, xe)


