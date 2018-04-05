import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman.mineral import Mineral
from util import BurnManTest
from burnman.tools import *

from burnman.material import material_property
from burnman.tools import copy_documentation


class test_decorators__material_property(BurnManTest):

    class MyCountingMaterial(burnman.Material):

        def __init__(self):
            burnman.Material.__init__(self)
            self.counter = 0

        @material_property
        def some_property(self):
            """my documentation"""
            self.counter += 1
            return 1.0

    def test_caching(self):
        m = self.MyCountingMaterial()
        self.assertEqual(m.counter, 0)
        self.assertEqual(m.some_property, 1.0)
        self.assertEqual(m.counter, 1)
        self.assertEqual(m.some_property, 1.0)
        self.assertEqual(m.counter, 1)

    def test_reset(self):
        m = self.MyCountingMaterial()
        self.assertEqual(len(m._cached), 0)
        self.assertEqual(m.some_property, 1.0)
        self.assertEqual(len(m._cached), 1)

        m.reset()

        self.assertEqual(len(m._cached), 0)
        self.assertEqual(m.some_property, 1.0)
        self.assertEqual(len(m._cached), 1)

        self.assertEqual(m.counter, 2)

    def test_doc(self):
        """make sure documentation is passed through with the new decorator"""

        self.assertEqual(
            self.MyCountingMaterial.some_property.__doc__, "my documentation")

        internal_energy_doc = burnman.Material.molar_internal_energy.__doc__
        self.assertTrue("internal energy of" in internal_energy_doc)

        self.assertEqual(self.MyCountingMaterial.molar_internal_energy.__doc__,
                         burnman.Material.molar_internal_energy.__doc__)

        pressure_doc = burnman.Material.pressure.__doc__
        self.assertTrue("current pressure" in pressure_doc)


class ClassA_for_copy_documentation(object):

    def my_func(self):
        """A.my_func doc"""
        pass

    def func_with_args(self, x, y):
        """doc for func_with_args"""
        return x + y


class ClassB_for_copy_documentation(object):

    @copy_documentation(ClassA_for_copy_documentation.my_func)
    def some_func(self):
        """B.some_func doc"""
        pass

    @copy_documentation(ClassA_for_copy_documentation.my_func)
    def some_func_without_doc(self):
        pass


class test_decorators_copy_documentation(BurnManTest):

    def test(self):
        self.assertEqual(
            ClassA_for_copy_documentation.my_func.__doc__, "A.my_func doc")
        self.assertEqual(ClassB_for_copy_documentation.some_func.__doc__,
                         "(copied from my_func):\nA.my_func doc\nB.some_func doc")
        self.assertEqual(
            ClassB_for_copy_documentation.some_func_without_doc.__doc__, "(copied from my_func):\nA.my_func doc")

    def test_with_args(self):
        """ early versions couldn't deal with functions with parameters"""
        class C(object):

            @copy_documentation(ClassA_for_copy_documentation.func_with_args)
            def another(self):
                return 1.0

        self.assertEqual(
            C.another.__doc__, "(copied from func_with_args):\ndoc for func_with_args")

        c = C()
        self.assertEqual(c.another(), 1.0)

    def test_with_args2(self):
        """ early versions couldn't deal with functions with parameters"""
        class C(object):

            @copy_documentation(ClassA_for_copy_documentation.func_with_args)
            def bla(self, z):
                "bla"
                return 1.0 + z

        self.assertEqual(
            C.bla.__doc__, "(copied from func_with_args):\ndoc for func_with_args\nbla")
        c = C()
        self.assertEqual(c.bla(1.0), 2.0)


class ClassC_for_copy_documentation(burnman.Material):

    def __init__(self):
        burnman.Material.__init__(self)

    @material_property
    @copy_documentation(ClassA_for_copy_documentation.my_func)
    def some_func(self):
        """C.some_func doc"""
        return 1.0

    @material_property
    def some_func2(self):
        return 2.0

    @material_property
    @copy_documentation(ClassA_for_copy_documentation.my_func)
    def some_func3(self):
        return 3.0


class test_two_decorators(BurnManTest):

    def check(self, C):
        self.assertEqual(C.some_func, 1.0)
        self.assertEqual(C.some_func2, 2.0)
        self.assertEqual(C.some_func3, 3.0)

    def test(self):
        c = ClassC_for_copy_documentation()

        self.check(c)
        self.check(c)
        c.reset()
        self.check(c)
        self.check(c)


if __name__ == '__main__':
    unittest.main()
