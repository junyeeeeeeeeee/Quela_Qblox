from test_folder.test1 import MyClass
from test_folder.test2 import function1
from test_folder import test3 as t3

def t4():
   print("t4 works!")

__all__ = [
   '__version__',
   'MyClass',
]