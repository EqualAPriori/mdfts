import mdfts.forcefield.forcefield as ff
import mdfts.utils.serial as serial

BeadTypeList = serial.serialize_list(ff.BeadType)
print('=== New class: ===')
print(BeadTypeList)
print(BeadTypeList.__name__)
btl = BeadTypeList()

print("\n=== Test 1: initialize from list ===")
print(">>> btl.from_dict([['A',1.0],['B',2.0]])")
btl.from_dict([['A',1.0],['B',2.0]])
print(btl)

print("\n=== Save params to dict ===")
print(">>> params = btl.to_dict()")
params = btl.to_dict()
print(params)


print("\n=== Test 2: initialize from params dict, should give same result as above ===")
print(">>> btl.from_dict(params)")
btl.from_dict(params)
print(btl)


print("\n=== Test 3: elements of input can themselves be a mix of lists and dictionaries ===")
print(">>> btl.from_dict([['A',1.0],{'name':'B','smear_length':2.0}])")
btl.from_dict([['A',1.0],{'name':'B','smear_length':2.0}])
print(btl)


print("\n=== Test 4: Using SerialzableTypedList ===")
BTL = serial.SerializableTypedList(ff.BeadType)
BTL.from_dict([['A',1.0],{'name':'B','smear_length':2.0}])
print(BTL)

