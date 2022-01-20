#import sys, os

#parent_dir = os.path.dirname(__file__)
#sys.path.insert(1, parent_dir)
from mdfts.utils.serial import Serializable, serialize


@serialize(["x"])
class f:
    def __init__(self, x):
        self.x = x


F = f(7)


# === Demonstrating NESTED serialization works, and with properties too! ===
@serialize(["temp", "d"])
class proptest1:
    _serial_vars = []
    # _serial_vars.append('temp')
    # _serial_vars.append('d')

    def __init__(self, t=0):
        self._temp = t

        # Add a serializable member/attribute
        self.d = Serializable()
        self.d.x = 5
        self.d._serial_vars.append("x")

    @property
    def temp(self):
        print("In getter")
        return self._temp

    @temp.setter
    def temp(self, value):
        print("In setter")
        self._temp = value


p1 = proptest1(37)
print(p1.to_dict())

