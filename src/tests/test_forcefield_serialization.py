import mdfts.forcefield.forcefield as ff


print('=== sample force field')
f = ff.ForceField()
f.add_bead_type( ff.BeadType('A',2.0,-1.0) )
f.add_bead_type( ff.BeadType('B',5.0,1.0) )
fdict = f.to_dict()
print(f)

print("\nserializes to:\n{}".format(fdict))


print("\n=== Make another forcefield from the above")
g = ff.ForceField()
g.from_dict(fdict)
gdict = g.to_dict()
print(g)
print("\nserializes to:\n{}".format(gdict))

