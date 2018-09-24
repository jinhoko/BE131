import sys

def addSwigInterfacePath(version=2):
    if version == 3:
        sys.path.insert(0, '../interfaces/Python3')
    else:
        sys.path.insert(0, '../interfaces/Python')

def getDataDirPath():
    return "../tests/data/"
