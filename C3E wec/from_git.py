import capytaine as cpt

wi = 25 # width of flap [m]
th = 1 # thickness of flap [m]
h = 19 # height of flap [m], draft = 16 m
x = 0
y = 0
draft = 16 # m
z = 0.5*h-draft # box center [m]
cog = -11.4 # center of gravity above still water line [m]
nw = 40 # number of panels along width (x)
nt = 2 # number of panels along thickness (y)
nh = 40 # number of panels along height (z)

body = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(wi, th, h),
                                                                             resolution=(nw, nt, nh),
                                                                             center=(x, y, z),
                                                                             name='rect'))
body.keep_immersed_part()
body.center_of_mass=(0,0,cog)
dofs = body.add_all_rigid_body_dofs()
body.inertia_matrix = body.compute_rigid_body_inertia()
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

# Check:
awl = wi * th  # m**2
kb = 0.5 * draft  # m
kg = draft+cog # m
bmt = th ** 2 / 12 / draft  # m
bml = wi ** 2 / 12 / draft  # m
gmt = kb + bmt - kg  # m
gml = kb + bml - kg  # m

g = 9.81 # m/s**2
rho = 1000 # kg/m**3
disp = wi * th * draft # m**3
m = disp * rho # kg
c33 = rho * g * awl # N/m
c44 = m * g * gmt # Nm/rad
c55 = m * g * gml # Nm/rad

stiffnes_matrix = body.hydrostatic_stiffness.values

print('C11: {:.2f} vs {:.2f}, ratio {:.2f}'.format(stiffnes_matrix[0,0],c33,stiffnes_matrix[0,0]/c33))
print('C33: {:.2f} vs {:.2f}, ratio {:.2f}'.format(stiffnes_matrix[2,2],c33,stiffnes_matrix[2,2]/c33))
print('C44: {:.2f} vs {:.2f}, ratio {:.2f}'.format(stiffnes_matrix[3,3],c44,stiffnes_matrix[3,3]/c44))
print('C55: {:.2f} vs {:.2f}, ratio {:.2f}'.format(stiffnes_matrix[4,4],c55,stiffnes_matrix[4,4]/c55))

print('stiffness matrix',stiffnes_matrix)

#body.show()