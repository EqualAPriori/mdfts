from __future__ import absolute_import, division, print_function

from mdfts.utils import serial


@serial.serialize([])
class _Potential(object):

    def __init__(self):
        pass
